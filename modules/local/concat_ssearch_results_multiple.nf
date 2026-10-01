process CONCAT_SSEARCH_RESULTS {
    tag { "${meta.id}" }
    label 'light'
    publishDir ({ "${params.outdir}/${meta.id}/ssearch" }, mode: 'copy', pattern: '*.out')

    input:
    tuple val(meta), path(chunk_results)

    output:
    tuple val(meta), path("${meta.id}.ssearch_fragment_realignment_results.out"), emit: ssearch_results

    script:
    """
    # Nextflow automatically expands the list of chunk paths into a space-separated string
    cat ${chunk_results} > ${meta.id}.ssearch_fragment_realignment_results.out
    """
}