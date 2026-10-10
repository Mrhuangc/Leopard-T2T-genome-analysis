process INTERPROSCAN_CANDIDATE_SCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../interproscan_functional_smoke/environment.yml"

    input:
    tuple val(meta), path(proteins)
    path smoke_gate
    val interproscan_dir
    val expected_version
    val expected_ipr
    val applications
    val disable_precalc

    output:
    tuple val(meta), path("${meta.id}.interproscan_candidates.tsv"), path("${meta.id}.interproscan_candidates.version.txt"), emit: result
    tuple val(meta), path("${meta.id}.interproscan_candidates.log"), emit: log
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    def disablePrecalc = disable_precalc ? '-dp' : ''
    """
    fail_preflight() {
        echo "FATAL: formal InterProScan preflight: \$*" >&2
        exit 1
    }

    [[ -s "$proteins" ]] || fail_preflight "candidate protein FASTA is empty"
    grep -Fqx 'status=pass' "$smoke_gate" || \
        fail_preflight "functional smoke gate is absent or did not pass"
    grep -Fqx "expected_interpro_accession=$expected_ipr" "$smoke_gate" || \
        fail_preflight "functional smoke gate used a different InterPro accession"
    [[ ",$applications," == *",PrositeProfiles,"* ]] || \
        fail_preflight "PrositeProfiles is required for $expected_ipr"
    [[ -x "${interproscan_dir}/interproscan.sh" ]] || \
        fail_preflight "not executable: ${interproscan_dir}/interproscan.sh"
    [[ -f "${interproscan_dir}/interproscan.properties" ]] || \
        fail_preflight "not found: ${interproscan_dir}/interproscan.properties"

    "${interproscan_dir}/interproscan.sh" --version \
        > "${prefix}.interproscan_candidates.version.txt" 2>&1
    grep -Fq "$expected_version" "${prefix}.interproscan_candidates.version.txt" || \
        fail_preflight "expected InterProScan version $expected_version"

    find_pftool() {
        local candidate resolved
        for candidate in "\$@"; do
            resolved=\$(command -v "\$candidate" 2>/dev/null || true)
            if [[ -n "\$resolved" && -x "\$resolved" ]]; then
                printf '%s\\n' "\$resolved"
                return 0
            fi
        done
        return 1
    }

    pfscan_path=\$(find_pftool pfscanV3 pfscanv3 pfscan || true)
    pfsearch_path=\$(find_pftool pfsearchV3 pfsearchv3 pfsearch || true)
    [[ -n "\$pfscan_path" ]] || fail_preflight "pftools pfscan not found in PATH"
    [[ -n "\$pfsearch_path" ]] || fail_preflight "pftools pfsearch not found in PATH"
    for pftool_path in "\$pfscan_path" "\$pfsearch_path"; do
        missing_libraries=\$(ldd "\$pftool_path" 2>&1 | awk '/not found/' || true)
        [[ -z "\$missing_libraries" ]] || \
            fail_preflight "missing shared libraries for \$pftool_path: \$missing_libraries"
    done

    mkdir -p "\$PWD/.interproscan-5"
    private_properties="\$PWD/.interproscan-5/interproscan.properties"
    cp "${interproscan_dir}/interproscan.properties" "\$private_properties"
    printf '\\n%s\\n%s\\n' \
        "binary.prosite.pfscanv3.path=\$pfscan_path" \
        "binary.prosite.pfsearchv3.path=\$pfsearch_path" \
        >> "\$private_properties"
    export _JAVA_OPTIONS="-Duser.home=\$PWD"

    mkdir -p interproscan_tmp
    scan_status=0
    "${interproscan_dir}/interproscan.sh" \
        -i "$proteins" \
        -appl "$applications" \
        -f TSV \
        -b "${prefix}.interproscan_candidates" \
        -cpu "$task.cpus" \
        -T "\$PWD/interproscan_tmp" \
        $disablePrecalc \
        > "${prefix}.interproscan_candidates.log" 2>&1 || scan_status=\$?

    if [[ "\$scan_status" -ne 0 ]]; then
        echo "FATAL: formal InterProScan exited with status \$scan_status" >&2
        echo '----- InterProScan log (last 200 lines) -----' >&2
        tail -n 200 "${prefix}.interproscan_candidates.log" >&2 || true
        exit "\$scan_status"
    fi
    [[ -f "${prefix}.interproscan_candidates.tsv" ]] || \
        fail_preflight "InterProScan returned success without a TSV output"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: $expected_version
        pftools: 3.2.13
        applications: $applications
        expected_interpro_accession: $expected_ipr
        java: \$(java -version 2>&1 | awk 'NR == 1 {gsub(/"/, "", \$3); print \$3}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    candidate_id=\$(awk '/^>/{sub(/^>/, ""); split(\$0, fields, /[[:space:]]+/); print fields[1]; exit}' "$proteins")
    echo -e "\${candidate_id}\\tmd5\\t300\\tProSiteProfiles\\tPS50262\\tGPCR profile\\t1\\t300\\t1e-20\\tT\\t2026-08-12\\t$expected_ipr\\tGPCR, rhodopsin-like, 7TM" > "${prefix}.interproscan_candidates.tsv"
    echo 'InterProScan version $expected_version' > "${prefix}.interproscan_candidates.version.txt"
    echo 'stub' > "${prefix}.interproscan_candidates.log"
    echo '"${task.process}": {interproscan: "$expected_version", pftools: "3.2.13"}' > versions.yml
    """
}
