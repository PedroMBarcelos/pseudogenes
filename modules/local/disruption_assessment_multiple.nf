process DISRUPTION_ASSESSMENT {
    tag { "${meta.id}" }
    label 'medium'
    cpus 4

    publishDir({ "${params.outdir}/${meta.id}/pseudogenes" }, mode: 'copy', overwrite: true)

    input:
    tuple val(meta), path(initial_jsonl), path(genome_fasta)
    path reference_fasta
    path assessment_script
    path tfasty_executable
    val flank_window
    val alignment_timeout

    output:
    tuple val(meta), path("${meta.id}.tfasty_v3_assessed.jsonl"), emit: pseudogenes
    tuple val(meta), path("${meta.id}.tfasty_v3_assessed.audit.tsv"), emit: audit
    tuple val(meta), path("${meta.id}.tfasty_v3_assessed.failed.jsonl"), emit: failures

    script:
    """
    python3 "${assessment_script}" \\
        "${initial_jsonl}" \\
        "${meta.id}.tfasty_v3_assessed.jsonl" \\
        --reference "${reference_fasta}" \\
        --genome "${genome_fasta}" \\
        --tfasty "${tfasty_executable}" \\
        --flank-window ${flank_window} \\
        --timeout ${alignment_timeout} \\
        --workers ${task.cpus}
    """
}
