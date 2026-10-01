process MAKE_MIN_EVALUE_FILE {
    label 'light'
    publishDir "${params.outdir}/metrics", mode: 'copy', pattern: 'min_evalue.txt'

    input:
    val min_evalue_value

    output:
    path 'min_evalue.txt', emit: min_evalue

    script:
    """
    printf "%s\n" "${min_evalue_value}" > min_evalue.txt
    """
}
