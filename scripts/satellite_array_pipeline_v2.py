#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
satellite_array_pipeline.py

针对高分歧、circular monomer 边界不稳定的着丝粒 satellite（例如猫科 SAT113）进行
array-specific 分析。

核心思路
--------
1. 不用单一全局 consensus 先过滤掉大部分 copy；保留全部 TRASH2 candidate。
2. 按 chromosome + TRASH2 array_id 分组。
3. 每个 array 选择局部 seed，第一次统一方向和 circular phase。
4. 每个 array 独立 MAFFT，生成第一轮 local consensus。
5. 再用 local consensus 对原始 copy 做第二轮方向/phase 标准化。
6. 再次 MAFFT，得到最终 array-specific MSA 和 consensus。
7. 输出：
   - 每个 copy 相对本阵列 consensus 的 substitution divergence
   - deletion-aware divergence
   - 阵列内 pairwise divergence（全部 copy 与 full-length-only 分开统计）
   - divergence 沿 TRASH2 block 位置的变化
   - 合并相邻 block 后，沿大型 SAT-rich region 的位置梯度
   - 先全局统一方向和 circular phase，再计算 array consensus 间距离

输入 TRASH2 repeats_with_seq CSV（无表头，11 列）：
chrom,array_id,start,end,strand,score,evalue,length,name,family_eds,sequence

推荐运行：
python satellite_array_pipeline.py \
  -i PPA.final.genome.fna_repeats_with_seq.SAT113.csv \
  -o PPA.SAT113.array_pipeline \
  --prefix PPA.SAT113 \
  --min-array-copies 5 \
  --threads 20

依赖：
- Python 3.8+
- MAFFT（必须在 PATH 中，或通过 --mafft 指定）

说明：
- 默认假定 repeats_with_seq CSV 的 sequence 列已经按 TRASH2 注释方向输出，因此
  --sequence-orientation 默认为 as_is。
- 如果你的 CSV sequence 实际是基因组正链序列，请使用：
  --sequence-orientation reverse_negative
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


VALID = set("ACGT")
ALN_VALID = set("ACGT-")
RC_TABLE = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")
CSV_COLUMNS = [
    "chrom", "array_id", "start", "end", "strand", "score", "evalue",
    "length", "name", "family_eds", "sequence",
]


@dataclass
class CopyRecord:
    copy_id: str
    chrom: str
    array_id: str
    array_key: str
    start: int
    end: int
    strand: str
    score: float
    evalue: float
    length: int
    name: str
    family_eds: str
    sequence: str
    row_number: int


@dataclass
class PhaseResult:
    copy_id: str
    orientation: str
    rotation: int
    identity: float
    matches: int
    compared: int
    votes: int
    sequence: str
    status: str


def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()


def clean_seq(seq: str) -> str:
    seq = seq.upper().replace(" ", "").replace("\t", "").replace("-", "")
    return "".join(base if base in VALID else "N" for base in seq)


def safe_name(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text))
    return text.strip("_") or "array"


def fmt_float(value: Optional[float], digits: int = 6) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{value:.{digits}f}"


def percentile(values: Sequence[float], q: float) -> float:
    vals = sorted(v for v in values if not math.isnan(v))
    if not vals:
        return math.nan
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - pos) + vals[hi] * (pos - lo)


def summary_stats(values: Sequence[float]) -> Dict[str, float]:
    vals = [v for v in values if not math.isnan(v)]
    if not vals:
        return {k: math.nan for k in ["mean", "median", "q25", "q75", "min", "max"]}
    return {
        "mean": statistics.fmean(vals),
        "median": statistics.median(vals),
        "q25": percentile(vals, 0.25),
        "q75": percentile(vals, 0.75),
        "min": min(vals),
        "max": max(vals),
    }


def parse_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_trash_csv(path: str, sequence_orientation: str) -> List[CopyRecord]:
    records: List[CopyRecord] = []
    seen_ids = Counter()

    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        for row_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if len(row) < 11:
                raise ValueError(
                    f"{path} 第 {row_number} 行只有 {len(row)} 列；需要至少 11 列。\n{row}"
                )

            chrom = row[0].strip()
            array_id = row[1].strip()
            start = int(row[2])
            end = int(row[3])
            strand = row[4].strip()
            score = parse_float(row[5])
            evalue = parse_float(row[6])
            length = int(row[7])
            name = row[8].strip()
            family_eds = row[9].strip()
            seq = clean_seq(row[10])

            if end <= start:
                raise ValueError(f"第 {row_number} 行坐标非法：start={start}, end={end}")
            if not seq:
                raise ValueError(f"第 {row_number} 行序列为空。")

            if sequence_orientation == "reverse_negative" and strand == "-":
                seq = revcomp(seq)

            array_key = f"{chrom}__{array_id}"
            base_id = f"{chrom}:{start}-{end}({strand})"
            seen_ids[base_id] += 1
            copy_id = base_id if seen_ids[base_id] == 1 else f"{base_id}__dup{seen_ids[base_id]}"

            records.append(
                CopyRecord(
                    copy_id=copy_id,
                    chrom=chrom,
                    array_id=array_id,
                    array_key=array_key,
                    start=start,
                    end=end,
                    strand=strand,
                    score=score,
                    evalue=evalue,
                    length=length,
                    name=name,
                    family_eds=family_eds,
                    sequence=seq,
                    row_number=row_number,
                )
            )

    if not records:
        raise ValueError(f"未从 {path} 读取到任何记录。")
    return records


def write_fasta(records: Iterable[Tuple[str, str]], path: Path, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    chunks: List[str] = []

    def flush() -> None:
        nonlocal header, chunks
        if header is not None:
            seq = "".join(chunks).strip().upper()
            if not seq:
                raise ValueError(f"空 FASTA 序列：{header} in {path}")
            records.append((header, seq))
        header = None
        chunks = []

    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].split()[0]
            else:
                if header is None:
                    raise ValueError(f"FASTA 格式错误：{path}")
                chunks.append(line)
        flush()
    return records


def choose_seed(records: List[CopyRecord]) -> CopyRecord:
    """
    选择第一轮局部 seed：
    1. 优先 modal length；
    2. 同长度时优先 TRASH score 高者；
    3. 再优先靠近 array 几何中心的 copy。

    第一轮 seed 只用于建立粗略 local consensus；第二轮会由 array-specific consensus 重校正。
    """
    lengths = [r.length for r in records]
    length_counts = Counter(lengths)
    modal_length = sorted(length_counts.items(), key=lambda x: (-x[1], abs(x[0] - statistics.median(lengths)), x[0]))[0][0]
    center = (min(r.start for r in records) + max(r.end for r in records)) / 2.0

    def key(r: CopyRecord):
        score = r.score if not math.isnan(r.score) else float("-inf")
        midpoint = (r.start + r.end) / 2.0
        return (
            abs(r.length - modal_length),
            -score,
            abs(midpoint - center),
            r.start,
        )

    return min(records, key=key)


def rotate_left(seq: str, offset: int) -> str:
    if not seq:
        return seq
    offset %= len(seq)
    return seq[offset:] + seq[:offset]


def build_kmer_index(ref: str, k: int, max_kmer_hits: int) -> Dict[str, List[int]]:
    ref = clean_seq(ref)
    L = len(ref)
    ref2 = ref + ref[:k - 1]
    tmp: Dict[str, List[int]] = defaultdict(list)
    for pos in range(L):
        kmer = ref2[pos:pos + k]
        if len(kmer) == k and set(kmer) <= VALID:
            tmp[kmer].append(pos)
    return {kmer: pos for kmer, pos in tmp.items() if len(pos) <= max_kmer_hits}


def vote_candidate_rotations(ref: str, seq: str, k: int, index: Dict[str, List[int]]) -> Counter:
    votes: Counter = Counter()
    n = len(seq)
    if n < k:
        return votes
    for i in range(0, n - k + 1):
        kmer = seq[i:i + k]
        if set(kmer) <= VALID and kmer in index:
            for ref_pos in index[kmer]:
                votes[(i - ref_pos) % n] += 1
    return votes


def get_candidate_offsets(
    votes: Counter,
    n: int,
    top_n: int,
    refine_window: int,
    fallback_max_rotations: int,
) -> List[int]:
    if n <= 0:
        return [0]
    if votes:
        candidates = set()
        for rot, _ in votes.most_common(top_n):
            for delta in range(-refine_window, refine_window + 1):
                candidates.add((rot + delta) % n)
        return sorted(candidates)

    # 没有 k-mer vote 时，不直接丢弃。
    # 对短 monomer 全部 rotation 穷举；较长序列则均匀抽样 rotation。
    if n <= fallback_max_rotations:
        return list(range(n))
    step = max(1, n // fallback_max_rotations)
    return list(range(0, n, step))[:fallback_max_rotations]


def ungapped_identity_to_ref(ref: str, seq: str, rotation: int) -> Tuple[float, int, int]:
    rotated = rotate_left(seq, rotation)
    m = min(len(ref), len(rotated))
    matches = 0
    compared = 0
    for a, b in zip(ref[:m], rotated[:m]):
        if a in VALID and b in VALID:
            compared += 1
            if a == b:
                matches += 1
    if compared == 0:
        return 0.0, 0, 0
    return matches / compared, matches, compared


def phase_one_sequence(
    ref: str,
    seq: str,
    k: int,
    max_kmer_hits: int,
    top_candidates: int,
    refine_window: int,
    min_identity: float,
    fallback_max_rotations: int,
) -> PhaseResult:
    ref = clean_seq(ref)
    seq = clean_seq(seq)
    index = build_kmer_index(ref, k, max_kmer_hits)
    best: Optional[dict] = None

    for orientation, oriented in [("forward", seq), ("reverse_complement", revcomp(seq))]:
        votes = vote_candidate_rotations(ref, oriented, k, index)
        offsets = get_candidate_offsets(
            votes, len(oriented), top_candidates, refine_window, fallback_max_rotations
        )
        for rot in offsets:
            ident, matches, compared = ungapped_identity_to_ref(ref, oriented, rot)
            result = {
                "orientation": orientation,
                "rotation": rot,
                "identity": ident,
                "matches": matches,
                "compared": compared,
                "votes": votes.get(rot, 0),
                "sequence": rotate_left(oriented, rot),
            }
            if best is None:
                best = result
            else:
                key = (result["identity"], result["votes"], result["matches"], -result["rotation"])
                best_key = (best["identity"], best["votes"], best["matches"], -best["rotation"])
                if key > best_key:
                    best = result

    if best is None:
        raise RuntimeError("无法完成 phase normalization。")

    return PhaseResult(
        copy_id="",
        orientation=best["orientation"],
        rotation=best["rotation"],
        identity=best["identity"],
        matches=best["matches"],
        compared=best["compared"],
        votes=best["votes"],
        sequence=best["sequence"],
        status="PASS" if best["identity"] >= min_identity else "LOW_CONFIDENCE",
    )


def phase_array(
    records: List[CopyRecord],
    ref_seq: str,
    out_fasta: Path,
    report_writer: csv.writer,
    stage: str,
    args,
) -> None:
    output: List[Tuple[str, str]] = []
    for rec in records:
        result = phase_one_sequence(
            ref=ref_seq,
            seq=rec.sequence,
            k=args.phase_k,
            max_kmer_hits=args.max_kmer_hits,
            top_candidates=args.top_candidates,
            refine_window=args.refine_window,
            min_identity=args.phase_min_identity,
            fallback_max_rotations=args.fallback_max_rotations,
        )
        result.copy_id = rec.copy_id
        output.append((rec.copy_id, result.sequence))
        report_writer.writerow([
            stage, rec.array_key, rec.copy_id, len(rec.sequence), len(ref_seq),
            result.orientation, result.rotation,
            fmt_float(result.identity * 100, 6), result.matches, result.compared,
            result.votes, result.status,
        ])
    write_fasta(output, out_fasta)


def run_mafft(
    mafft: str, infile: Path, outfile: Path, threads: int, quiet: bool = False,
    resume: bool = False,
) -> None:
    outfile.parent.mkdir(parents=True, exist_ok=True)
    if resume and outfile.exists() and outfile.stat().st_size > 0:
        print(f"      [resume] reuse MAFFT: {outfile}", file=sys.stderr)
        return
    cmd = [mafft, "--auto", "--thread", str(threads), str(infile)]
    stderr = subprocess.DEVNULL if quiet else None
    with open(outfile, "w", encoding="utf-8") as out:
        proc = subprocess.run(cmd, stdout=out, stderr=stderr, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"MAFFT 运行失败（exit={proc.returncode}）：{' '.join(cmd)}")


def build_consensus_from_alignment(
    aln_records: List[Tuple[str, str]],
    min_column_occupancy: float,
    min_consensus_support: float,
) -> Tuple[str, str, dict]:
    """
    返回：
    - ungapped consensus（低 occupancy 列删除）
    - aligned consensus（与 MSA 等长；删除列用 '-'）
    - stats
    """
    if not aln_records:
        raise ValueError("空 alignment。")
    lengths = {len(seq) for _, seq in aln_records}
    if len(lengths) != 1:
        raise ValueError("alignment 中序列长度不一致。")

    n = len(aln_records)
    L = lengths.pop()
    kept_chars: List[str] = []
    aligned_chars: List[str] = []
    kept_columns = 0
    low_occupancy_columns = 0
    n_columns_with_N = 0

    for idx in range(L):
        counts = Counter(seq[idx].upper() for _, seq in aln_records if seq[idx].upper() in VALID)
        canonical = sum(counts.values())
        occupancy = canonical / n
        if canonical == 0 or occupancy < min_column_occupancy:
            aligned_chars.append("-")
            low_occupancy_columns += 1
            continue
        base, count = sorted(counts.items(), key=lambda x: (-x[1], x[0]))[0]
        support = count / canonical
        if support < min_consensus_support:
            base = "N"
            n_columns_with_N += 1
        aligned_chars.append(base)
        kept_chars.append(base)
        kept_columns += 1

    consensus = "".join(kept_chars)
    aligned_consensus = "".join(aligned_chars)
    stats = {
        "n_sequences": n,
        "alignment_length": L,
        "consensus_length": len(consensus),
        "kept_columns": kept_columns,
        "low_occupancy_columns": low_occupancy_columns,
        "consensus_N_columns": n_columns_with_N,
    }
    return consensus, aligned_consensus, stats


def copy_vs_consensus_metrics(seq: str, aligned_consensus: str) -> dict:
    if len(seq) != len(aligned_consensus):
        raise ValueError("copy 与 aligned consensus 长度不一致。")

    sub_comp = sub_mismatch = 0
    del_comp = del_diff = 0
    consensus_bases = 0
    covered_consensus_bases = 0

    for a, c in zip(seq.upper(), aligned_consensus.upper()):
        if c in VALID:
            consensus_bases += 1
            if a in VALID:
                covered_consensus_bases += 1
                sub_comp += 1
                if a != c:
                    sub_mismatch += 1
            if a in VALID or a == "-":
                del_comp += 1
                if a != c:
                    del_diff += 1

    sub_div = 100.0 * sub_mismatch / sub_comp if sub_comp else math.nan
    del_div = 100.0 * del_diff / del_comp if del_comp else math.nan
    coverage = covered_consensus_bases / consensus_bases if consensus_bases else math.nan

    return {
        "substitution_comparable_sites": sub_comp,
        "substitution_mismatches": sub_mismatch,
        "substitution_divergence_pct": sub_div,
        "deletion_aware_comparable_sites": del_comp,
        "deletion_aware_differences": del_diff,
        "deletion_aware_divergence_pct": del_div,
        "consensus_coverage": coverage,
    }


def pair_divergence(seq1: str, seq2: str) -> Tuple[float, float, int, int]:
    if len(seq1) != len(seq2):
        raise ValueError("pairwise alignment 长度不一致。")

    sub_comp = sub_diff = 0
    indel_comp = indel_diff = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a in VALID and b in VALID:
            sub_comp += 1
            if a != b:
                sub_diff += 1
        if a in ALN_VALID and b in ALN_VALID and not (a == "-" and b == "-"):
            indel_comp += 1
            if a != b:
                indel_diff += 1

    sub = 100.0 * sub_diff / sub_comp if sub_comp else math.nan
    indel = 100.0 * indel_diff / indel_comp if indel_comp else math.nan
    return sub, indel, sub_comp, indel_comp


def choose_pairs(n: int, max_pairs: int, rng: random.Random) -> List[Tuple[int, int]]:
    total = n * (n - 1) // 2
    if total <= max_pairs:
        return list(itertools.combinations(range(n), 2))

    chosen = set()
    while len(chosen) < max_pairs:
        i = rng.randrange(n)
        j = rng.randrange(n - 1)
        if j >= i:
            j += 1
        if i > j:
            i, j = j, i
        chosen.add((i, j))
    return sorted(chosen)


def quality_tier(
    raw_length: int,
    consensus_length: int,
    consensus_coverage: float,
    phase_status: str,
    min_ratio: float,
    max_ratio: float,
    min_coverage: float,
) -> str:
    ratio = raw_length / consensus_length if consensus_length > 0 else math.nan
    if phase_status == "LOW_CONFIDENCE":
        return "phase_low_confidence"
    if not math.isnan(ratio) and min_ratio <= ratio <= max_ratio and consensus_coverage >= min_coverage:
        return "full_length"
    return "fragment_or_extended"


def inter_consensus_divergence(aln_path: Path, out_long: Path, out_matrix: Path) -> None:
    records = read_fasta(aln_path)
    names = [h for h, _ in records]
    seqs = [s for _, s in records]

    matrix: Dict[Tuple[str, str], float] = {}
    with open(out_long, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "array1", "array2", "substitution_divergence_pct",
            "indel_aware_divergence_pct", "substitution_comparable_sites",
            "indel_aware_comparable_sites",
        ])
        for i in range(len(records)):
            matrix[(names[i], names[i])] = 0.0
            for j in range(i + 1, len(records)):
                sub, indel, sub_comp, indel_comp = pair_divergence(seqs[i], seqs[j])
                matrix[(names[i], names[j])] = sub
                matrix[(names[j], names[i])] = sub
                writer.writerow([
                    names[i], names[j], fmt_float(sub, 6), fmt_float(indel, 6),
                    sub_comp, indel_comp,
                ])

    with open(out_matrix, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["array"] + names)
        for a in names:
            writer.writerow([a] + [fmt_float(matrix.get((a, b), math.nan), 6) for b in names])



def pairwise_divergence_summary(
    aln_records: List[Tuple[str, str]],
    max_pairs: int,
    rng: random.Random,
) -> Tuple[Dict[str, float], Dict[str, float], int]:
    """对一组已在同一 MSA 中的序列计算抽样/全量 pairwise divergence。"""
    pair_indices = choose_pairs(len(aln_records), max_pairs, rng)
    pair_sub: List[float] = []
    pair_indel: List[float] = []
    for i, j in pair_indices:
        sub, indel, _, _ = pair_divergence(aln_records[i][1], aln_records[j][1])
        if not math.isnan(sub):
            pair_sub.append(sub)
        if not math.isnan(indel):
            pair_indel.append(indel)
    return summary_stats(pair_sub), summary_stats(pair_indel), len(pair_indices)


def build_merged_regions(
    arrays: Dict[str, List[CopyRecord]],
    max_gap: int,
) -> Tuple[Dict[str, str], Dict[str, dict]]:
    """
    将同一染色体上相距不超过 max_gap 的 TRASH2 blocks 合并为 SAT-rich region。
    返回 array_key -> region_key 映射和 region 元数据。
    """
    by_chrom: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for key, arr in arrays.items():
        by_chrom[arr[0].chrom].append((min(r.start for r in arr), max(r.end for r in arr), key))

    array_to_region: Dict[str, str] = {}
    regions: Dict[str, dict] = {}
    for chrom in sorted(by_chrom):
        blocks = sorted(by_chrom[chrom])
        groups: List[dict] = []
        current: Optional[dict] = None
        for start, end, key in blocks:
            if current is None or start - current["end"] > max_gap:
                if current is not None:
                    groups.append(current)
                current = {"start": start, "end": end, "array_keys": [key]}
            else:
                current["end"] = max(current["end"], end)
                current["array_keys"].append(key)
        if current is not None:
            groups.append(current)

        for idx, group in enumerate(groups, start=1):
            region_key = f"{chrom}__region{idx:03d}"
            region = {
                "region_key": region_key,
                "chrom": chrom,
                "region_start": group["start"],
                "region_end": group["end"],
                "region_span_bp": group["end"] - group["start"],
                "n_arrays": len(group["array_keys"]),
                "member_arrays": sorted(group["array_keys"]),
            }
            regions[region_key] = region
            for key in group["array_keys"]:
                array_to_region[key] = region_key
    return array_to_region, regions


def choose_array_consensus_anchor(
    consensuses: Dict[str, str],
    arrays: Dict[str, List[CopyRecord]],
    consensus_stats: Dict[str, dict],
) -> str:
    """选择长度典型、N 少且 copy 数较多的 array consensus 作为全局 phase anchor。"""
    if not consensuses:
        raise ValueError("没有 array consensus 可用于选择全局 anchor。")
    lengths = [len(seq) for seq in consensuses.values()]
    modal_length = Counter(lengths).most_common(1)[0][0]

    def rank(key: str):
        seq = consensuses[key]
        n_count = seq.count("N")
        return (
            abs(len(seq) - modal_length),
            n_count / max(1, len(seq)),
            -len(arrays.get(key, [])),
            key,
        )

    return min(consensuses, key=rank)


def phase_consensus_collection(
    consensuses: Dict[str, str],
    ref: str,
    stage: str,
    args,
    report_writer: csv.writer,
) -> Dict[str, str]:
    normalized: Dict[str, str] = {}
    for key in sorted(consensuses):
        seq = consensuses[key]
        result = phase_one_sequence(
            ref=ref,
            seq=seq,
            k=args.phase_k,
            max_kmer_hits=args.max_kmer_hits,
            top_candidates=args.top_candidates,
            refine_window=args.refine_window,
            min_identity=args.inter_array_phase_min_identity,
            fallback_max_rotations=args.fallback_max_rotations,
        )
        normalized[key] = result.sequence
        report_writer.writerow([
            stage, key, len(seq), len(ref), result.orientation, result.rotation,
            fmt_float(result.identity * 100, 6), result.matches, result.compared,
            result.votes, result.status,
        ])
    return normalized


def write_position_bins(
    rows: List[dict],
    group_field: str,
    frac_field: str,
    output_path: Path,
    position_bins: int,
) -> None:
    by_group_bin: Dict[Tuple[str, int], List[dict]] = defaultdict(list)
    for row in rows:
        frac = min(max(float(row[frac_field]), 0.0), 1.0)
        b = min(position_bins - 1, int(frac * position_bins))
        by_group_bin[(str(row[group_field]), b)].append(row)

    with open(output_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            group_field, "chrom", "bin_index", "bin_start_frac", "bin_end_frac",
            "n_copies", "n_full_length", "full_length_fraction",
            "median_substitution_divergence_all_pct",
            "median_substitution_divergence_full_length_pct",
            "median_deletion_aware_divergence_all_pct",
            "median_deletion_aware_divergence_full_length_pct",
        ])
        for (group, b), group_rows in sorted(by_group_bin.items()):
            all_sub = [float(r["substitution_divergence_pct"]) for r in group_rows
                       if not math.isnan(float(r["substitution_divergence_pct"]))]
            all_del = [float(r["deletion_aware_divergence_pct"]) for r in group_rows
                       if not math.isnan(float(r["deletion_aware_divergence_pct"]))]
            full_rows = [r for r in group_rows if r["quality_tier"] == "full_length"]
            full_sub = [float(r["substitution_divergence_pct"]) for r in full_rows
                        if not math.isnan(float(r["substitution_divergence_pct"]))]
            full_del = [float(r["deletion_aware_divergence_pct"]) for r in full_rows
                        if not math.isnan(float(r["deletion_aware_divergence_pct"]))]
            writer.writerow([
                group, group_rows[0]["chrom"], b + 1,
                fmt_float(b / position_bins, 6), fmt_float((b + 1) / position_bins, 6),
                len(group_rows), len(full_rows), fmt_float(len(full_rows) / len(group_rows), 6),
                fmt_float(statistics.median(all_sub) if all_sub else math.nan, 6),
                fmt_float(statistics.median(full_sub) if full_sub else math.nan, 6),
                fmt_float(statistics.median(all_del) if all_del else math.nan, 6),
                fmt_float(statistics.median(full_del) if full_del else math.nan, 6),
            ])


def write_metadata(
    records: List[CopyRecord],
    arrays: Dict[str, List[CopyRecord]],
    seeds: Dict[str, CopyRecord],
    meta_dir: Path,
    min_array_copies: int,
) -> None:
    meta_dir.mkdir(parents=True, exist_ok=True)

    with open(meta_dir / "master_candidates.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "copy_id", "chrom", "array_id", "array_key", "region_key",
        "region_start", "region_end", "region_span_bp", "start", "end", "strand",
            "score", "evalue", "reported_length", "sequence_length", "trash_name",
            "family_eds", "source_row",
        ])
        for r in records:
            writer.writerow([
                r.copy_id, r.chrom, r.array_id, r.array_key, r.start, r.end, r.strand,
                fmt_float(r.score, 6), fmt_float(r.evalue, 6), r.length, len(r.sequence),
                r.name, r.family_eds, r.row_number,
            ])

    with open(meta_dir / "array_manifest.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "array_key", "array_filename", "chrom", "array_id", "n_copies",
            "array_start", "array_end", "array_span_bp", "median_copy_length",
            "modal_copy_length", "min_copy_length", "max_copy_length", "trash_names",
            "seed_copy_id", "seed_length", "analyzed",
        ])
        for key in sorted(arrays):
            arr = arrays[key]
            lengths = [r.length for r in arr]
            mode_len = Counter(lengths).most_common(1)[0][0]
            seed = seeds[key]
            writer.writerow([
                key, safe_name(key), arr[0].chrom, arr[0].array_id, len(arr),
                min(r.start for r in arr), max(r.end for r in arr),
                max(r.end for r in arr) - min(r.start for r in arr),
                fmt_float(statistics.median(lengths), 3), mode_len, min(lengths), max(lengths),
                ",".join(sorted(set(r.name for r in arr))), seed.copy_id, seed.length,
                "yes" if len(arr) >= min_array_copies else "no",
            ])


def main() -> None:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Array-specific pipeline for highly divergent centromeric satellite monomers.",
    )
    parser.add_argument("-i", "--input-csv", required=True, help="TRASH2 repeats_with_seq CSV（11 列，无表头）")
    parser.add_argument("-o", "--output-dir", required=True, help="输出目录")
    parser.add_argument("--prefix", default="SAT", help="结果前缀")
    parser.add_argument("--min-array-copies", type=int, default=5, help="进入 array-specific MSA 的最少 copy 数")
    parser.add_argument("--threads", type=int, default=8, help="MAFFT 线程数")
    parser.add_argument("--mafft", default="mafft", help="MAFFT 可执行文件路径")
    parser.add_argument(
        "--sequence-orientation", choices=["as_is", "reverse_negative"], default="as_is",
        help="CSV sequence 是否直接使用；若 sequence 是基因组正链，选 reverse_negative",
    )

    # phase normalization
    parser.add_argument("--phase-k", type=int, default=9, help="local phase normalization 使用的 k-mer 长度")
    parser.add_argument("--max-kmer-hits", type=int, default=20, help="reference 中出现过多次的 k-mer 被忽略")
    parser.add_argument("--top-candidates", type=int, default=20, help="每个方向保留的 rotation 候选数")
    parser.add_argument("--refine-window", type=int, default=5, help="候选 rotation 左右微调窗口")
    parser.add_argument("--phase-min-identity", type=float, default=0.35, help="仅用于标记 LOW_CONFIDENCE，不删除 copy")
    parser.add_argument("--fallback-max-rotations", type=int, default=300, help="无 k-mer vote 时最多检查的 rotation 数")

    # consensus
    parser.add_argument("--min-column-occupancy", type=float, default=0.50, help="建立 local consensus 时保留 alignment column 的最小 occupancy")
    parser.add_argument("--min-consensus-support", type=float, default=0.50, help="major allele 低于该比例时 consensus 写 N")

    # full-length label, not hard filtering
    parser.add_argument("--full-length-min-ratio", type=float, default=0.80, help="raw length/local consensus length 最小值")
    parser.add_argument("--full-length-max-ratio", type=float, default=1.20, help="raw length/local consensus length 最大值")
    parser.add_argument("--full-length-min-coverage", type=float, default=0.80, help="copy 覆盖 local consensus 的最小比例")

    # pairwise / position
    parser.add_argument("--max-pairwise-pairs", type=int, default=10000, help="每个 array 最多计算的 copy pair 数")
    parser.add_argument("--random-seed", type=int, default=20260715, help="pairwise 抽样随机种子")
    parser.add_argument("--position-bins", type=int, default=20, help="沿 array/region 位置分箱数")
    parser.add_argument(
        "--merge-block-gap", type=int, default=10000,
        help="同一染色体相邻 TRASH2 blocks 间距不超过该值时合并为 SAT-rich region；设为 0 仅合并重叠/相邻 block",
    )
    parser.add_argument(
        "--inter-array-phase-min-identity", type=float, default=0.30,
        help="array consensus 全局方向/循环相位校正的 LOW_CONFIDENCE 标记阈值（不删除）",
    )

    parser.add_argument("--prepare-only", action="store_true", help="只解析 CSV 并输出 metadata/raw FASTA，不运行 MAFFT")
    parser.add_argument(
        "--resume", action="store_true",
        help="若已有非空 MAFFT 输出则复用，可在原输出目录上快速重算统计与修正后的阵列间距离",
    )
    parser.add_argument("--quiet-mafft", action="store_true", help="隐藏 MAFFT stderr")
    args = parser.parse_args()

    if args.min_array_copies < 2:
        parser.error("--min-array-copies 至少为 2。")
    if args.phase_k < 5:
        parser.error("--phase-k 建议至少为 5。")
    if not (0 <= args.phase_min_identity <= 1):
        parser.error("--phase-min-identity 必须在 0–1。")
    if not (0 < args.min_column_occupancy <= 1):
        parser.error("--min-column-occupancy 必须在 (0,1]。")
    if not (0 < args.min_consensus_support <= 1):
        parser.error("--min-consensus-support 必须在 (0,1]。")
    if args.merge_block_gap < 0:
        parser.error("--merge-block-gap 不能小于 0。")
    if not (0 <= args.inter_array_phase_min_identity <= 1):
        parser.error("--inter-array-phase-min-identity 必须在 0–1。")

    outdir = Path(args.output_dir)
    meta_dir = outdir / "00_metadata"
    raw_dir = outdir / "01_raw_arrays"
    pass1_phase_dir = outdir / "02_pass1_phase"
    pass1_aln_dir = outdir / "03_pass1_mafft"
    local_cons_dir = outdir / "04_local_consensus_pass1"
    final_phase_dir = outdir / "05_final_phase"
    final_aln_dir = outdir / "06_final_mafft"
    result_dir = outdir / "07_results"

    for d in [meta_dir, raw_dir, pass1_phase_dir, pass1_aln_dir, local_cons_dir, final_phase_dir, final_aln_dir, result_dir]:
        d.mkdir(parents=True, exist_ok=True)

    print(f"[1/9] Reading candidate CSV: {args.input_csv}", file=sys.stderr)
    records = read_trash_csv(args.input_csv, args.sequence_orientation)
    arrays: Dict[str, List[CopyRecord]] = defaultdict(list)
    rec_by_id: Dict[str, CopyRecord] = {}
    for r in records:
        arrays[r.array_key].append(r)
        rec_by_id[r.copy_id] = r
    for key in arrays:
        arrays[key].sort(key=lambda x: (x.start, x.end, x.copy_id))

    seeds = {key: choose_seed(arr) for key, arr in arrays.items()}
    write_metadata(records, arrays, seeds, meta_dir, args.min_array_copies)

    analyzed_arrays = {k: v for k, v in arrays.items() if len(v) >= args.min_array_copies}
    excluded_arrays = {k: v for k, v in arrays.items() if len(v) < args.min_array_copies}
    array_to_region, merged_regions = build_merged_regions(analyzed_arrays, args.merge_block_gap)

    with open(meta_dir / "merged_region_manifest.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "region_key", "chrom", "region_start", "region_end", "region_span_bp",
            "n_arrays", "member_arrays", "merge_block_gap",
        ])
        for region_key in sorted(merged_regions):
            region = merged_regions[region_key]
            writer.writerow([
                region_key, region["chrom"], region["region_start"], region["region_end"],
                region["region_span_bp"], region["n_arrays"],
                ",".join(region["member_arrays"]), args.merge_block_gap,
            ])

    for key, arr in arrays.items():
        write_fasta([(r.copy_id, r.sequence) for r in arr], raw_dir / f"{safe_name(key)}.raw.fa")
    write_fasta([(r.copy_id, r.sequence) for r in records], raw_dir / f"{args.prefix}.all_candidates.raw.fa")

    print(
        f"      candidates={len(records)}; arrays={len(arrays)}; "
        f"analyzed_arrays={len(analyzed_arrays)}; excluded_arrays={len(excluded_arrays)}; "
        f"merged_regions={len(merged_regions)}",
        file=sys.stderr,
    )

    if args.prepare_only:
        print("[done] --prepare-only requested.", file=sys.stderr)
        return

    mafft_path = shutil.which(args.mafft) if os.path.sep not in args.mafft else args.mafft
    if not mafft_path or not Path(mafft_path).exists():
        raise SystemExit(
            f"找不到 MAFFT：{args.mafft}\n"
            "请先安装 MAFFT，或使用 --mafft /path/to/mafft。"
        )

    # phase report files
    pass1_phase_report = result_dir / f"{args.prefix}.pass1.phase_report.tsv"
    final_phase_report = result_dir / f"{args.prefix}.final.phase_report.tsv"
    phase_header = [
        "stage", "array_key", "copy_id", "input_length", "reference_length",
        "orientation", "rotation_left_bp", "ungapped_identity_to_reference_pct",
        "matches", "compared_sites", "kmer_votes", "phase_status",
    ]

    print("[2/9] First-pass array-local phase normalization", file=sys.stderr)
    with open(pass1_phase_report, "w", encoding="utf-8", newline="") as rep_fh:
        rep = csv.writer(rep_fh, delimiter="\t")
        rep.writerow(phase_header)
        for idx, key in enumerate(sorted(analyzed_arrays), start=1):
            arr = analyzed_arrays[key]
            seed = seeds[key]
            phase_array(
                arr, seed.sequence,
                pass1_phase_dir / f"{safe_name(key)}.phase1.fa",
                rep, "pass1_seed", args,
            )
            if idx % 10 == 0 or idx == len(analyzed_arrays):
                print(f"      phased {idx}/{len(analyzed_arrays)} arrays", file=sys.stderr)

    print("[3/9] First-pass MAFFT per array", file=sys.stderr)
    for idx, key in enumerate(sorted(analyzed_arrays), start=1):
        run_mafft(
            mafft=str(mafft_path),
            infile=pass1_phase_dir / f"{safe_name(key)}.phase1.fa",
            outfile=pass1_aln_dir / f"{safe_name(key)}.pass1.aln.fa",
            threads=args.threads,
            quiet=args.quiet_mafft, resume=args.resume,
        )
        if idx % 10 == 0 or idx == len(analyzed_arrays):
            print(f"      aligned {idx}/{len(analyzed_arrays)} arrays", file=sys.stderr)

    print("[4/9] Build pass1 array-specific local consensuses", file=sys.stderr)
    local_consensus: Dict[str, str] = {}
    with open(result_dir / f"{args.prefix}.pass1.consensus_stats.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "array_key", "n_sequences", "alignment_length", "consensus_length",
            "kept_columns", "low_occupancy_columns", "consensus_N_columns",
        ])
        for key in sorted(analyzed_arrays):
            aln = read_fasta(pass1_aln_dir / f"{safe_name(key)}.pass1.aln.fa")
            cons, _, stats = build_consensus_from_alignment(
                aln, args.min_column_occupancy, args.min_consensus_support
            )
            if not cons:
                raise RuntimeError(f"{key} 第一轮 local consensus 为空。")
            local_consensus[key] = cons
            write_fasta([(key, cons)], local_cons_dir / f"{safe_name(key)}.local_consensus.fa")
            writer.writerow([
                key, stats["n_sequences"], stats["alignment_length"], stats["consensus_length"],
                stats["kept_columns"], stats["low_occupancy_columns"], stats["consensus_N_columns"],
            ])
    write_fasta(sorted(local_consensus.items()), local_cons_dir / f"{args.prefix}.all_array_consensus.pass1.fa")

    print("[5/9] Second-pass phase normalization against local consensus", file=sys.stderr)
    final_phase_status: Dict[str, str] = {}
    with open(final_phase_report, "w", encoding="utf-8", newline="") as rep_fh:
        rep = csv.writer(rep_fh, delimiter="\t")
        rep.writerow(phase_header)
        for idx, key in enumerate(sorted(analyzed_arrays), start=1):
            arr = analyzed_arrays[key]
            ref = local_consensus[key]
            # phase_array 写报告，同时下面重新获取 status 供 quality tier 使用
            output: List[Tuple[str, str]] = []
            for rec in arr:
                result = phase_one_sequence(
                    ref=ref,
                    seq=rec.sequence,
                    k=args.phase_k,
                    max_kmer_hits=args.max_kmer_hits,
                    top_candidates=args.top_candidates,
                    refine_window=args.refine_window,
                    min_identity=args.phase_min_identity,
                    fallback_max_rotations=args.fallback_max_rotations,
                )
                final_phase_status[rec.copy_id] = result.status
                output.append((rec.copy_id, result.sequence))
                rep.writerow([
                    "final_local_consensus", key, rec.copy_id, len(rec.sequence), len(ref),
                    result.orientation, result.rotation,
                    fmt_float(result.identity * 100, 6), result.matches, result.compared,
                    result.votes, result.status,
                ])
            write_fasta(output, final_phase_dir / f"{safe_name(key)}.phase2.fa")
            if idx % 10 == 0 or idx == len(analyzed_arrays):
                print(f"      re-phased {idx}/{len(analyzed_arrays)} arrays", file=sys.stderr)

    print("[6/9] Final MAFFT per array", file=sys.stderr)
    for idx, key in enumerate(sorted(analyzed_arrays), start=1):
        run_mafft(
            mafft=str(mafft_path),
            infile=final_phase_dir / f"{safe_name(key)}.phase2.fa",
            outfile=final_aln_dir / f"{safe_name(key)}.final.aln.fa",
            threads=args.threads,
            quiet=args.quiet_mafft, resume=args.resume,
        )
        if idx % 10 == 0 or idx == len(analyzed_arrays):
            print(f"      aligned {idx}/{len(analyzed_arrays)} arrays", file=sys.stderr)

    print("[7/9] Calculate copy-to-consensus, intra-array and positional divergence", file=sys.stderr)
    copy_rows: List[dict] = []
    array_summary_rows: List[dict] = []
    position_rows: List[dict] = []
    final_array_consensus: Dict[str, str] = {}
    final_array_consensus_stats: Dict[str, dict] = {}
    rng = random.Random(args.random_seed)

    final_cons_stats_path = result_dir / f"{args.prefix}.final.consensus_stats.tsv"
    with open(final_cons_stats_path, "w", encoding="utf-8", newline="") as cons_stats_fh:
        cons_writer = csv.writer(cons_stats_fh, delimiter="\t")
        cons_writer.writerow([
            "array_key", "n_sequences", "alignment_length", "consensus_length",
            "kept_columns", "low_occupancy_columns", "consensus_N_columns",
        ])

        for idx, key in enumerate(sorted(analyzed_arrays), start=1):
            arr = analyzed_arrays[key]
            aln_records = read_fasta(final_aln_dir / f"{safe_name(key)}.final.aln.fa")
            aln_dict = dict(aln_records)
            cons, aligned_cons, cons_stats = build_consensus_from_alignment(
                aln_records, args.min_column_occupancy, args.min_consensus_support
            )
            final_array_consensus[key] = cons
            final_array_consensus_stats[key] = cons_stats
            cons_writer.writerow([
                key, cons_stats["n_sequences"], cons_stats["alignment_length"],
                cons_stats["consensus_length"], cons_stats["kept_columns"],
                cons_stats["low_occupancy_columns"], cons_stats["consensus_N_columns"],
            ])

            arr_start = min(r.start for r in arr)
            arr_end = max(r.end for r in arr)
            arr_span = max(1, arr_end - arr_start)
            sorted_arr = sorted(arr, key=lambda r: (r.start, r.end, r.copy_id))
            n_arr = len(sorted_arr)

            sub_divs_all: List[float] = []
            del_divs_all: List[float] = []
            sub_divs_full: List[float] = []
            del_divs_full: List[float] = []
            tier_by_id: Dict[str, str] = {}

            for rank, rec in enumerate(sorted_arr):
                seq = aln_dict.get(rec.copy_id)
                if seq is None:
                    raise KeyError(f"{key} alignment 中找不到 copy：{rec.copy_id}")
                m = copy_vs_consensus_metrics(seq, aligned_cons)
                midpoint = (rec.start + rec.end) / 2.0
                genomic_pos_frac = (midpoint - arr_start) / arr_span
                rank_frac = rank / (n_arr - 1) if n_arr > 1 else 0.5
                edge_distance_rank_frac = 2.0 * min(rank_frac, 1.0 - rank_frac)
                ratio = rec.length / len(cons) if len(cons) > 0 else math.nan
                pstatus = final_phase_status.get(rec.copy_id, "NA")
                tier = quality_tier(
                    rec.length, len(cons), m["consensus_coverage"], pstatus,
                    args.full_length_min_ratio, args.full_length_max_ratio,
                    args.full_length_min_coverage,
                )

                region_key = array_to_region[key]
                region = merged_regions[region_key]
                row = {
                    "copy_id": rec.copy_id,
                    "chrom": rec.chrom,
                    "array_id": rec.array_id,
                    "array_key": key,
                    "region_key": region_key,
                    "region_start": region["region_start"],
                    "region_end": region["region_end"],
                    "region_span_bp": region["region_span_bp"],
                    "start": rec.start,
                    "end": rec.end,
                    "strand": rec.strand,
                    "trash_name": rec.name,
                    "raw_length": rec.length,
                    "local_consensus_length": len(cons),
                    "length_ratio_to_local_consensus": ratio,
                    "phase_status": pstatus,
                    "quality_tier": tier,
                    **m,
                    "array_rank": rank + 1,
                    "array_rank_frac": rank_frac,
                    "array_genomic_pos_frac": genomic_pos_frac,
                    "edge_to_center_rank_score": edge_distance_rank_frac,
                }
                tier_by_id[rec.copy_id] = tier
                copy_rows.append(row)
                position_rows.append(row)
                sub_divs_all.append(m["substitution_divergence_pct"])
                del_divs_all.append(m["deletion_aware_divergence_pct"])
                if tier == "full_length":
                    sub_divs_full.append(m["substitution_divergence_pct"])
                    del_divs_full.append(m["deletion_aware_divergence_pct"])

            # intra-array pairwise divergence：全部 copy 与 full-length-only 分开统计。
            pair_sub_all_stats, pair_indel_all_stats, n_pairs_all = pairwise_divergence_summary(
                aln_records, args.max_pairwise_pairs, rng
            )
            full_aln_records = [
                (copy_id, seq) for copy_id, seq in aln_records
                if tier_by_id.get(copy_id) == "full_length"
            ]
            pair_sub_full_stats, pair_indel_full_stats, n_pairs_full = pairwise_divergence_summary(
                full_aln_records, args.max_pairwise_pairs, rng
            )

            sub_all_stats = summary_stats(sub_divs_all)
            del_all_stats = summary_stats(del_divs_all)
            sub_full_stats = summary_stats(sub_divs_full)
            del_full_stats = summary_stats(del_divs_full)

            array_summary_rows.append({
                "array_key": key,
                "chrom": arr[0].chrom,
                "array_id": arr[0].array_id,
                "n_candidates": len(arr),
                "n_full_length": len(sub_divs_full),
                "full_length_fraction": len(sub_divs_full) / len(arr),
                "array_start": arr_start,
                "array_end": arr_end,
                "array_span_bp": arr_end - arr_start,
                "local_consensus_length": len(cons),
                "copy_to_consensus_sub_div_mean_all_pct": sub_all_stats["mean"],
                "copy_to_consensus_sub_div_median_all_pct": sub_all_stats["median"],
                "copy_to_consensus_sub_div_q25_all_pct": sub_all_stats["q25"],
                "copy_to_consensus_sub_div_q75_all_pct": sub_all_stats["q75"],
                "copy_to_consensus_del_div_median_all_pct": del_all_stats["median"],
                "copy_to_consensus_sub_div_median_full_length_pct": sub_full_stats["median"],
                "copy_to_consensus_del_div_median_full_length_pct": del_full_stats["median"],
                # 兼容旧列名：以下旧列均代表 all copies。
                "pairwise_pairs_evaluated": n_pairs_all,
                "intra_array_pairwise_sub_div_mean_pct": pair_sub_all_stats["mean"],
                "intra_array_pairwise_sub_div_median_pct": pair_sub_all_stats["median"],
                "intra_array_pairwise_sub_div_q25_pct": pair_sub_all_stats["q25"],
                "intra_array_pairwise_sub_div_q75_pct": pair_sub_all_stats["q75"],
                "intra_array_pairwise_indel_div_median_pct": pair_indel_all_stats["median"],
                # 新增显式 all/full-length-only 指标。
                "pairwise_pairs_evaluated_all": n_pairs_all,
                "intra_array_pairwise_sub_div_mean_all_pct": pair_sub_all_stats["mean"],
                "intra_array_pairwise_sub_div_median_all_pct": pair_sub_all_stats["median"],
                "intra_array_pairwise_sub_div_q25_all_pct": pair_sub_all_stats["q25"],
                "intra_array_pairwise_sub_div_q75_all_pct": pair_sub_all_stats["q75"],
                "intra_array_pairwise_indel_div_median_all_pct": pair_indel_all_stats["median"],
                "pairwise_pairs_evaluated_full_length": n_pairs_full,
                "intra_array_pairwise_sub_div_mean_full_length_pct": pair_sub_full_stats["mean"],
                "intra_array_pairwise_sub_div_median_full_length_pct": pair_sub_full_stats["median"],
                "intra_array_pairwise_sub_div_q25_full_length_pct": pair_sub_full_stats["q25"],
                "intra_array_pairwise_sub_div_q75_full_length_pct": pair_sub_full_stats["q75"],
                "intra_array_pairwise_indel_div_median_full_length_pct": pair_indel_full_stats["median"],
            })

            if idx % 10 == 0 or idx == len(analyzed_arrays):
                print(f"      analyzed divergence {idx}/{len(analyzed_arrays)} arrays", file=sys.stderr)

    # 为合并后的 SAT-rich region 重新计算全局位置 rank。
    rows_by_region: Dict[str, List[dict]] = defaultdict(list)
    for row in copy_rows:
        rows_by_region[row["region_key"]].append(row)
    for region_key, region_rows in rows_by_region.items():
        region = merged_regions[region_key]
        ordered = sorted(region_rows, key=lambda r: (int(r["start"]), int(r["end"]), r["copy_id"]))
        n_region = len(ordered)
        span = max(1, region["region_span_bp"])
        for rank, row in enumerate(ordered):
            midpoint = (int(row["start"]) + int(row["end"])) / 2.0
            rank_frac = rank / (n_region - 1) if n_region > 1 else 0.5
            row["region_rank"] = rank + 1
            row["region_rank_frac"] = rank_frac
            row["region_genomic_pos_frac"] = (midpoint - region["region_start"]) / span
            row["region_edge_to_center_rank_score"] = 2.0 * min(rank_frac, 1.0 - rank_frac)

    # region-level summary
    region_summary_path = result_dir / f"{args.prefix}.region_divergence_summary.tsv"
    with open(region_summary_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "region_key", "chrom", "region_start", "region_end", "region_span_bp",
            "n_arrays", "member_arrays", "n_candidates", "n_full_length",
            "full_length_fraction", "median_substitution_divergence_all_pct",
            "median_substitution_divergence_full_length_pct",
            "median_deletion_aware_divergence_all_pct",
            "median_deletion_aware_divergence_full_length_pct",
        ])
        for region_key in sorted(rows_by_region):
            rows = rows_by_region[region_key]
            region = merged_regions[region_key]
            full = [r for r in rows if r["quality_tier"] == "full_length"]
            sub_all = [float(r["substitution_divergence_pct"]) for r in rows]
            sub_full = [float(r["substitution_divergence_pct"]) for r in full]
            del_all = [float(r["deletion_aware_divergence_pct"]) for r in rows]
            del_full = [float(r["deletion_aware_divergence_pct"]) for r in full]
            writer.writerow([
                region_key, region["chrom"], region["region_start"], region["region_end"],
                region["region_span_bp"], region["n_arrays"], ",".join(region["member_arrays"]),
                len(rows), len(full), fmt_float(len(full) / len(rows), 6),
                fmt_float(statistics.median(sub_all) if sub_all else math.nan, 6),
                fmt_float(statistics.median(sub_full) if sub_full else math.nan, 6),
                fmt_float(statistics.median(del_all) if del_all else math.nan, 6),
                fmt_float(statistics.median(del_full) if del_full else math.nan, 6),
            ])

    # write copy table
    copy_columns = [
        "copy_id", "chrom", "array_id", "array_key", "region_key",
        "region_start", "region_end", "region_span_bp", "start", "end", "strand",
        "trash_name", "raw_length", "local_consensus_length",
        "length_ratio_to_local_consensus", "phase_status", "quality_tier",
        "substitution_comparable_sites", "substitution_mismatches",
        "substitution_divergence_pct", "deletion_aware_comparable_sites",
        "deletion_aware_differences", "deletion_aware_divergence_pct",
        "consensus_coverage", "array_rank", "array_rank_frac",
        "array_genomic_pos_frac", "edge_to_center_rank_score",
        "region_rank", "region_rank_frac", "region_genomic_pos_frac",
        "region_edge_to_center_rank_score",
    ]
    copy_div_path = result_dir / f"{args.prefix}.copy_divergence.tsv"
    with open(copy_div_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=copy_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in copy_rows:
            out = row.copy()
            for col in [
                "length_ratio_to_local_consensus", "substitution_divergence_pct",
                "deletion_aware_divergence_pct", "consensus_coverage", "array_rank_frac",
                "array_genomic_pos_frac", "edge_to_center_rank_score",
                "region_rank_frac", "region_genomic_pos_frac",
                "region_edge_to_center_rank_score",
            ]:
                out[col] = fmt_float(out[col], 6)
            writer.writerow(out)

    array_summary_path = result_dir / f"{args.prefix}.array_divergence_summary.tsv"
    array_summary_columns = list(array_summary_rows[0].keys()) if array_summary_rows else []
    with open(array_summary_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=array_summary_columns, delimiter="\t")
        writer.writeheader()
        for row in array_summary_rows:
            out = row.copy()
            for k, v in list(out.items()):
                if isinstance(v, float):
                    out[k] = fmt_float(v, 6)
            writer.writerow(out)

    # positional bins：分别输出 TRASH2 block 和合并 region，且区分 all/full-length。
    bins_path = result_dir / f"{args.prefix}.array_position_bins.tsv"
    region_bins_path = result_dir / f"{args.prefix}.region_position_bins.tsv"
    write_position_bins(
        position_rows, "array_key", "array_rank_frac", bins_path, args.position_bins
    )
    write_position_bins(
        position_rows, "region_key", "region_rank_frac", region_bins_path, args.position_bins
    )

    # final array consensus：先全局统一方向和 circular phase，再进行最终 MAFFT。
    array_cons_raw_fa = result_dir / f"{args.prefix}.array_consensus.raw.fa"
    array_cons_phase1_fa = result_dir / f"{args.prefix}.array_consensus.global_phase1.fa"
    array_cons_phase1_aln = result_dir / f"{args.prefix}.array_consensus.global_phase1.aln.fa"
    array_cons_fa = result_dir / f"{args.prefix}.array_consensus.fa"
    array_cons_aln = result_dir / f"{args.prefix}.array_consensus.aln.fa"
    array_cons_phase_report = result_dir / f"{args.prefix}.array_consensus.phase_report.tsv"
    global_cons_fa = result_dir / f"{args.prefix}.global_SAT_consensus.fa"
    write_fasta(sorted(final_array_consensus.items()), array_cons_raw_fa)

    print("[8/9] Global orientation/circular-phase normalization of array consensuses", file=sys.stderr)
    anchor_key = choose_array_consensus_anchor(
        final_array_consensus, analyzed_arrays, final_array_consensus_stats
    )
    anchor_seq = final_array_consensus[anchor_key]
    with open(array_cons_phase_report, "w", encoding="utf-8", newline="") as fh:
        reporter = csv.writer(fh, delimiter="\t")
        reporter.writerow([
            "stage", "array_key", "input_length", "reference_length", "orientation",
            "rotation_left_bp", "ungapped_identity_to_reference_pct", "matches",
            "compared_sites", "kmer_votes", "phase_status",
        ])
        phase1_cons = phase_consensus_collection(
            final_array_consensus, anchor_seq, "anchor_phase", args, reporter
        )
        write_fasta(sorted(phase1_cons.items()), array_cons_phase1_fa)
        run_mafft(
            mafft=str(mafft_path), infile=array_cons_phase1_fa, outfile=array_cons_phase1_aln,
            threads=args.threads, quiet=args.quiet_mafft, resume=False,
        )
        phase1_aln_records = read_fasta(array_cons_phase1_aln)
        global_cons1, _, _ = build_consensus_from_alignment(
            phase1_aln_records, args.min_column_occupancy, args.min_consensus_support
        )
        phase2_cons = phase_consensus_collection(
            final_array_consensus, global_cons1, "global_consensus_phase", args, reporter
        )

    # 保留旧文件名，但内容现在是已全局 phase-normalized 的 consensus。
    write_fasta(sorted(phase2_cons.items()), array_cons_fa)
    print("[9/9] Final inter-array alignment and divergence", file=sys.stderr)
    run_mafft(
        mafft=str(mafft_path), infile=array_cons_fa, outfile=array_cons_aln,
        threads=args.threads, quiet=args.quiet_mafft, resume=False,
    )
    final_global_cons, _, final_global_stats = build_consensus_from_alignment(
        read_fasta(array_cons_aln), args.min_column_occupancy, args.min_consensus_support
    )
    write_fasta([(f"{args.prefix}_global_consensus", final_global_cons)], global_cons_fa)
    inter_consensus_divergence(
        array_cons_aln,
        result_dir / f"{args.prefix}.inter_array_divergence.long.tsv",
        result_dir / f"{args.prefix}.inter_array_substitution_divergence.matrix.tsv",
    )

    # run summary
    full_count = sum(r["quality_tier"] == "full_length" for r in copy_rows)
    low_phase = sum(r["phase_status"] == "LOW_CONFIDENCE" for r in copy_rows)
    with open(result_dir / f"{args.prefix}.run_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["input_candidates", len(records)])
        writer.writerow(["input_arrays", len(arrays)])
        writer.writerow(["analyzed_arrays", len(analyzed_arrays)])
        writer.writerow(["excluded_small_arrays", len(excluded_arrays)])
        writer.writerow(["merged_regions", len(merged_regions)])
        writer.writerow(["merge_block_gap", args.merge_block_gap])
        writer.writerow(["resume_mode", "yes" if args.resume else "no"])
        writer.writerow(["analyzed_copies", len(copy_rows)])
        writer.writerow(["full_length_copies", full_count])
        writer.writerow(["full_length_fraction", fmt_float(full_count / len(copy_rows) if copy_rows else math.nan, 6)])
        writer.writerow(["final_phase_low_confidence_copies", low_phase])
        writer.writerow(["max_pairwise_pairs_per_array", args.max_pairwise_pairs])
        writer.writerow(["position_bins", args.position_bins])
        writer.writerow(["inter_array_phase_anchor", anchor_key])
        writer.writerow(["global_consensus_length", len(final_global_cons)])
        writer.writerow(["global_consensus_N_columns", final_global_stats["consensus_N_columns"]])

    print("[done] Pipeline completed.", file=sys.stderr)
    print(f"       Main results: {result_dir}", file=sys.stderr)
    print(f"       Copy divergence: {copy_div_path}", file=sys.stderr)
    print(f"       Array summary: {array_summary_path}", file=sys.stderr)
    print(f"       Array position bins: {bins_path}", file=sys.stderr)
    print(f"       Region position bins: {region_bins_path}", file=sys.stderr)
    print(f"       Corrected inter-array matrix: {result_dir / f'{args.prefix}.inter_array_substitution_divergence.matrix.tsv'}", file=sys.stderr)


if __name__ == "__main__":
    main()
