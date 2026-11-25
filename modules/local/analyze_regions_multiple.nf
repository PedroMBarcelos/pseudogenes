process ANALYZE_REGIONS {
    tag { "${meta.id}" }
    label 'medium'

    publishDir(
        { "${params.outdir}/${meta.id}/analysis" },
        mode: 'copy',
        pattern: '*.tsv'
    )

    input:
    tuple val(meta),
          path(ssearch_results),
          path(sorted_bed),
          path(merged_bed),
          path(filtered_tblastn)

    path min_evalue

    output:
    tuple val(meta),
          path("${meta.id}.pseudogene_analysis_report.tsv"),
          emit: final_report

    script:
    """
    python3 "${projectDir}/scripts/analyze_regions_v2.py" \
      --ssearch "${ssearch_results}" \
      --tblastn "${filtered_tblastn}" \
      --bed_hits "${sorted_bed}" \
      --merged_regions "${merged_bed}" \
      --output "${meta.id}.pseudogene_analysis_report.tsv" \
      --min_evalue "\$(cat "${min_evalue}")"
    """
}