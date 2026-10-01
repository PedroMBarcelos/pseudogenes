process CONCAT_SSEARCH_RESULTS {
    label 'light'
    publishDir "${params.outdir}/ssearch", mode: 'copy', pattern: 'ssearch_fragment_realignment_results.out'

    input:
    path chunk_results

    output:
    path 'ssearch_fragment_realignment_results.out', emit: ssearch_results

    script:
    """
    cat ${chunk_results.collect { "\"${it}\"" }.join(' ')} > ssearch_fragment_realignment_results.out
    """
}
