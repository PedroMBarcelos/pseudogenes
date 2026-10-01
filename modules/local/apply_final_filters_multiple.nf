process APPLY_FINAL_FILTERS {
    tag "$meta.id"
    label 'light'
    publishDir "${params.outdir}/${meta.id}/final", mode: 'copy', pattern: '*.jsonl'
    
    input:
    tuple val(meta), path(pseudogene_annotations)
    val min_coverage
    val final_min_identity

    output:
    tuple val(meta), path("${meta.id}.pseudogene_annotations_with_coverage.jsonl"), emit: filtered_annotations

    script:
    """
    python3 ${workflow.projectDir}/apply_final_filters.py \
        "${pseudogene_annotations}" \
        "${min_coverage}" \
        "${final_min_identity}" \
        ./

    # Rename the output to include the unique genome ID to prevent race conditions
    if [ -f "pseudogene_annotations_with_coverage.jsonl" ]; then
        mv pseudogene_annotations_with_coverage.jsonl ${meta.id}.pseudogene_annotations_with_coverage.jsonl
    fi
    """
}