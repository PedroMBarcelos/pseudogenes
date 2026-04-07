process EXTRACT_MIN_EVALUE {
    label 'light'
    publishDir "${params.outdir}/metrics", mode: 'copy'

    input:
    path blastp_out

    output:
    path 'min_evalue.txt', emit: min_evalue

    script:
    """
    grep -v '^#' "${blastp_out}" | sort -k11,11g | head -n 1 | awk '{print \$11}' > min_evalue.txt
    """
}
