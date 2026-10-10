#define _GNU_SOURCE
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Deterministic global protein aligner for the FCA ChrY curation workflow.
 * Input:  pair_id<TAB>query_protein<TAB>target_protein
 * Output: pair_id, score, alignment length, matches, aligned residue pairs,
 *         gaps, identity, query coverage, target coverage, status.
 *
 * This deliberately uses a simple linear scoring scheme.  The aligner is
 * used only after a protein k-mer prefilter has selected close candidate
 * homologues, so the reported identity and two-sided coverage are the
 * relevant curation statistics.
 */

static int align_pair(const char *q, const char *t, FILE *out,
                      const char *pair_id) {
    const int MATCH = 2, MISMATCH = -1, GAP = -2;
    size_t n = strlen(q), m = strlen(t), cols = m + 1;
    if (n == 0 || m == 0) {
        fprintf(out, "%s\tNA\t0\t0\t0\t0\t0\t0\t0\tempty_sequence\n",
                pair_id);
        return 0;
    }
    if ((n + 1) > SIZE_MAX / (m + 1)) return 2;
    size_t cells = (n + 1) * (m + 1);
    if (cells > 30000000ULL) {
        fprintf(out, "%s\tNA\t0\t0\t0\t0\t0\t0\t0\ttoo_large\n",
                pair_id);
        return 0;
    }
    int32_t *score = (int32_t *)malloc(cells * sizeof(int32_t));
    uint8_t *trace = (uint8_t *)malloc(cells);
    if (!score || !trace) {
        free(score); free(trace);
        return 3;
    }
    score[0] = 0; trace[0] = 0;
    for (size_t i = 1; i <= n; ++i) {
        score[i * cols] = (int32_t)(i * GAP); trace[i * cols] = 1;
    }
    for (size_t j = 1; j <= m; ++j) {
        score[j] = (int32_t)(j * GAP); trace[j] = 2;
    }
    for (size_t i = 1; i <= n; ++i) {
        char a = (char)toupper((unsigned char)q[i - 1]);
        for (size_t j = 1; j <= m; ++j) {
            char b = (char)toupper((unsigned char)t[j - 1]);
            int32_t diag = score[(i - 1) * cols + j - 1] +
                           ((a == b && a != 'X') ? MATCH : MISMATCH);
            int32_t up = score[(i - 1) * cols + j] + GAP;
            int32_t left = score[i * cols + j - 1] + GAP;
            int32_t best = diag; uint8_t direction = 0;
            if (up > best) { best = up; direction = 1; }
            if (left > best) { best = left; direction = 2; }
            score[i * cols + j] = best;
            trace[i * cols + j] = direction;
        }
    }
    size_t i = n, j = m, alen = 0, matches = 0, both = 0, gaps = 0;
    while (i > 0 || j > 0) {
        uint8_t direction = trace[i * cols + j];
        if (i > 0 && j > 0 && direction == 0) {
            char a = (char)toupper((unsigned char)q[i - 1]);
            char b = (char)toupper((unsigned char)t[j - 1]);
            if (a == b && a != 'X') ++matches;
            ++both; --i; --j;
        } else if (i > 0 && (j == 0 || direction == 1)) {
            ++gaps; --i;
        } else {
            ++gaps; --j;
        }
        ++alen;
    }
    double identity = both ? (double)matches / (double)both : 0.0;
    double qcov = (double)both / (double)n;
    double tcov = (double)both / (double)m;
    fprintf(out, "%s\t%d\t%zu\t%zu\t%zu\t%zu\t%.6f\t%.6f\t%.6f\tok\n",
            pair_id, score[n * cols + m], alen, matches, both, gaps,
            identity, qcov, tcov);
    free(score); free(trace);
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s pairs.tsv results.tsv\n", argv[0]);
        return 2;
    }
    FILE *in = fopen(argv[1], "r");
    FILE *out = fopen(argv[2], "w");
    if (!in || !out) {
        fprintf(stderr, "Cannot open input or output file.\n");
        return 3;
    }
    fprintf(out, "pair_id\tscore\talignment_length\tmatches\taligned_residue_pairs\tgaps\tidentity\tquery_coverage\ttarget_coverage\tstatus\n");
    char *line = NULL; size_t cap = 0; ssize_t len;
    while ((len = getline(&line, &cap, in)) >= 0) {
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
            line[--len] = '\0';
        if (len == 0 || line[0] == '#') continue;
        char *p1 = strchr(line, '\t');
        if (!p1) continue;
        *p1++ = '\0';
        char *p2 = strchr(p1, '\t');
        if (!p2) continue;
        *p2++ = '\0';
        int rc = align_pair(p1, p2, out, line);
        if (rc != 0) {
            fprintf(stderr, "Alignment failed for %s (code %d).\n", line, rc);
            free(line); fclose(in); fclose(out); return 4;
        }
    }
    free(line); fclose(in); fclose(out);
    return 0;
}
