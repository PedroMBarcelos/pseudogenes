process PSEUDOGENE_DETECTION {
    tag { "${meta.id}" }
    label 'medium'
    publishDir({ "${params.outdir}/${meta.id}/pseudogenes" }, mode: 'copy', pattern: '*.jsonl')
    
    input:
    tuple val(meta), path(region_annotations), path(genome_fasta)
    
    output:
    tuple val(meta), path("${meta.id}.initial_pseudogene_annotations.jsonl"), emit: pseudogenes

    script:
    """
    python3 ${workflow.projectDir}/pseudogene_detection.py \
        "${region_annotations}" \
        "${genome_fasta}" \
        "${meta.id}.initial_pseudogene_annotations.jsonl"
    """
}