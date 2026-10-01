process FILTER_TBLASTN {
    tag { "${meta.id}" }
    label 'light'
    publishDir({ "${params.outdir}/${meta.id}/intermediate" }, mode: 'copy', pattern: '*.filtered.out')

    input:
    tuple val(meta), path(tblastn_out)  // <-- Must be a tuple here!
    path min_evalue

    output:
    tuple val(meta), path("${meta.id}.tblastn.filtered.out"), emit: filtered_tblastn

    script:
    """
    threshold=\$(cat "${min_evalue}")
    awk -v threshold="\${threshold}" '\$11 <= threshold' "${tblastn_out}" > ${meta.id}.tblastn.filtered.out
    """
}