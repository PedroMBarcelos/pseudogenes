process WRITE_FINAL_FILES {
    tag "$meta.id"
    label 'light'
    publishDir({ "${params.outdir}/${meta.id}/final" }, mode: 'copy', pattern: '*.{tsv,gff3,txt}')
    
    input:
    tuple val(meta), path(pseudogene_annotations)
    
    output:
    tuple val(meta), path("${meta.id}.pseudogene_annotations_final.tsv"), emit: final_tsv
    tuple val(meta), path("${meta.id}.pseudogene_annotations.gff3"), emit: final_gff
    tuple val(meta), path("${meta.id}.coverage_statistics.tsv"), emit: coverage_stats
    tuple val(meta), path("${meta.id}.summary_report.txt"), emit: summary_report

    script:
    """
    python3 ${workflow.projectDir}/write_final_files.py "${pseudogene_annotations}" ./
    
    # Safely rename outputs to include the unique genome ID and prevent race conditions
    [ -f "pseudogene_annotations_final.tsv" ] && mv pseudogene_annotations_final.tsv ${meta.id}.pseudogene_annotations_final.tsv
    [ -f "pseudogene_annotations.gff3" ] && mv pseudogene_annotations.gff3 ${meta.id}.pseudogene_annotations.gff3
    [ -f "coverage_statistics.tsv" ] && mv coverage_statistics.tsv ${meta.id}.coverage_statistics.tsv
    [ -f "summary_report.txt" ] && mv summary_report.txt ${meta.id}.summary_report.txt
    """
}