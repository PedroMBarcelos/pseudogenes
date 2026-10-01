process PREPARE_GENOME_FASTA {
    cache "lenient"
    tag "$meta.id"
    label 'light'
    publishDir "${params.outdir}/${meta.id}/intermediate", mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), path(genome_input), path(gff_file)

    output:
    tuple val(meta), path("${meta.id}.fasta"), path(gff_file), emit: genome_and_gff

    script:
    """
    if [[ "${genome_input}" == *.gz ]]; then
      gunzip -c "${genome_input}" > ${meta.id}.fasta
    else
      cp "${genome_input}" ${meta.id}.fasta
    fi
    """
}