process ANNOTATE_REGIONS_WITH_BEST_PROTEINS {
    tag { "${meta.id}" }
    label 'medium'

    publishDir(
        { "${params.outdir}/${meta.id}/analysis" },
        mode: 'copy',
        pattern: '*.region_annotations.jsonl'
    )

    input:
    tuple val(meta),
          path(analysis_report),
          path(filtered_tblastn),
          path(ssearch_results)

    output:
    tuple val(meta),
          path("${meta.id}.region_annotations.jsonl"),
          emit: annotations

    script:
    """
    python3 "${projectDir}/scripts/annotate_regions_with_best_proteins.py" \
        --regions "${analysis_report}" \
        --tblastn "${filtered_tblastn}" \
        --ssearch "${ssearch_results}" \
        --output "${meta.id}.region_annotations.jsonl"
    """
}