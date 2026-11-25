process WRITE_FINAL_FILES {
    tag "${meta.id}"
    label 'light'

    publishDir "${params.outdir}/${meta.id}/final",
        mode: 'copy'

    input:
    tuple val(meta), path(annotations_jsonl)

    output:
    tuple val(meta),
          path("${meta.id}.candidate_regions.tsv"),
          emit: final_tsv

    script:
    """
    python3 "${projectDir}/scripts/write_final_files.py" \\
        "${annotations_jsonl}" \\
        "${meta.id}.candidate_regions.tsv"
    """
}