process ANNOTATE_REGIONS_WITH_BEST_PROTEINS {
    tag { "${meta.id}" }
    label 'light'
    publishDir({ "${params.outdir}/${meta.id}/annotations" }, mode: 'copy', pattern: '*.jsonl')
    
    input:
    tuple val(meta), path(regions_path), path(hsps_path)
    
    output:
    tuple val(meta), path("${meta.id}.region_annotations.jsonl"), emit: annotations

    script:
    """
    python3 ${workflow.projectDir}/annotate_regions_with_best_proteins.py \
        "${regions_path}" \
        "${hsps_path}" \
        "${meta.id}.region_annotations.jsonl"
    """
}