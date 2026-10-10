process INTERPROSCAN_FUNCTIONAL_SMOKE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(proteins)
    val interproscan_dir
    val expected_version
    val applications
    val disable_precalc

    output:
    tuple val(meta), path("${meta.id}.interproscan_smoke.tsv"), path("${meta.id}.interproscan.version.txt"), emit: result
    tuple val(meta), path("${meta.id}.interproscan_smoke.log"), emit: log
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.id
    def disablePrecalc = disable_precalc ? '-dp' : ''
    """
    fail_preflight() {
        echo "FATAL: InterProScan preflight: \$*" >&2
        exit 1
    }

    [[ -x "${interproscan_dir}/interproscan.sh" ]] || \
        fail_preflight "not executable: ${interproscan_dir}/interproscan.sh"
    [[ -f "${interproscan_dir}/interproscan.properties" ]] || \
        fail_preflight "not found: ${interproscan_dir}/interproscan.properties"
    "${interproscan_dir}/interproscan.sh" --version > "${prefix}.interproscan.version.txt" 2>&1
    grep -Fq "$expected_version" "${prefix}.interproscan.version.txt" || \
        fail_preflight "expected version '$expected_version'; observed: \$(tr '\\n' ' ' < '${prefix}.interproscan.version.txt')"

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

    # Bioconda builds have used both versioned and unversioned command names.
    pfscan_path=\$(find_pftool pfscanV3 pfscanv3 pfscan || true)
    pfsearch_path=\$(find_pftool pfsearchV3 pfsearchv3 pfsearch || true)
    [[ -n "\$pfscan_path" ]] || \
        fail_preflight "pftools pfscan executable not found in PATH=\$PATH"
    [[ -n "\$pfsearch_path" ]] || \
        fail_preflight "pftools pfsearch executable not found in PATH=\$PATH"

    for pftool_path in "\$pfscan_path" "\$pfsearch_path"; do
        missing_libraries=\$(ldd "\$pftool_path" 2>&1 | awk '/not found/' || true)
        [[ -z "\$missing_libraries" ]] || \
            fail_preflight "missing shared libraries for \$pftool_path: \$missing_libraries"
    done

    echo "InterProScan preflight: pfscan=\$pfscan_path" >&2
    echo "InterProScan preflight: pfsearch=\$pfsearch_path" >&2

    mkdir -p "\$PWD/.interproscan-5"
    private_properties="\$PWD/.interproscan-5/interproscan.properties"
    cp "${interproscan_dir}/interproscan.properties" "\$private_properties"
    # java.util.Properties uses the last value for a duplicated key. Appending
    # task-local overrides avoids modifying the shared InterProScan install and
    # is robust to whitespace/layout changes in release property files.
    printf '\\n%s\\n%s\\n' \
        "binary.prosite.pfscanv3.path=\$pfscan_path" \
        "binary.prosite.pfsearchv3.path=\$pfsearch_path" \
        >> "\$private_properties"
    tail -n 2 "\$private_properties" >&2
    export _JAVA_OPTIONS="-Duser.home=\$PWD"

    mkdir -p interproscan_tmp
    scan_status=0
    "${interproscan_dir}/interproscan.sh" \
        -i "$proteins" \
        -appl "$applications" \
        -f TSV \
        -b "${prefix}.interproscan_smoke" \
        -cpu "$task.cpus" \
        -T "\$PWD/interproscan_tmp" \
        $disablePrecalc \
        > "${prefix}.interproscan_smoke.log" 2>&1 || scan_status=\$?

    if [[ "\$scan_status" -ne 0 ]]; then
        echo "FATAL: InterProScan exited with status \$scan_status" >&2
        echo "----- InterProScan log (last 200 lines) -----" >&2
        tail -n 200 "${prefix}.interproscan_smoke.log" >&2 || true
        exit "\$scan_status"
    fi

    [[ -s "${prefix}.interproscan_smoke.tsv" ]] || \
        fail_preflight "InterProScan returned success but the smoke-test TSV is missing or empty; inspect ${prefix}.interproscan_smoke.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: $expected_version
        pftools: 3.2.13
        applications: $applications
        java: \$(java -version 2>&1 | awk 'NR == 1 {gsub(/"/, "", \$3); print \$3}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    echo -e 'known_OR_smoke\tmd5\t5\tProSiteProfiles\tPS50262\tG-protein coupled receptors family 1 profile\t1\t5\t1e-10\tT\t2026-08-12\tIPR017452\tGPCR, rhodopsin-like, 7TM' > "${prefix}.interproscan_smoke.tsv"
    echo 'InterProScan version $expected_version' > "${prefix}.interproscan.version.txt"
    echo 'stub' > "${prefix}.interproscan_smoke.log"
    echo '"${task.process}": {interproscan: "$expected_version", pftools: "3.2.13", applications: "$applications", java: "11"}' > versions.yml
    """
}
