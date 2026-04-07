process ANALYZE_REGIONS {
    label 'medium'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'pseudogene_analysis_report.tsv'

    input:
    path ssearch_results
    path sorted_bed
    path merged_bed
    path min_evalue

    output:
    path 'pseudogene_analysis_report.tsv', emit: final_report

    script:
    """
    python3 "${projectDir}/analyze_regions.py" \
      --ssearch "${ssearch_results}" \
      --bed_hits "${sorted_bed}" \
      --merged_regions "${merged_bed}" \
      --output pseudogene_analysis_report.tsv \
      --min_evalue "\$(cat ${min_evalue})"
    """
}
