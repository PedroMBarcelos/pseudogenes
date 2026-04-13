process FILTER_TBLASTN {
    label 'light'
    publishDir "${params.outdir}/intermediate", mode: 'copy', pattern: 'genome.tblastn.filtered.out'

    input:
    path tblastn_out
    path min_evalue

    output:
    path 'genome.tblastn.filtered.out', emit: filtered_tblastn

    script:
    """
    threshold=\$(cat "${min_evalue}")
    awk -v threshold="\${threshold}" '\$11 <= threshold' "${tblastn_out}" > genome.tblastn.filtered.out
    """
}
