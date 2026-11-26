process WRITE_FINAL_FILES {
    tag "${meta.id}"
    label 'light'

    publishDir "${params.outdir}/${meta.id}/final",
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(annotations_jsonl)
    path exporter_script

    output:
    tuple val(meta),
          path("${meta.id}.candidate_regions.tsv"),
          emit: final_tsv

    script:
    """
    python3 "${exporter_script}" \\
        "${annotations_jsonl}" \\
        "${meta.id}.candidate_regions.tsv"
    """
}