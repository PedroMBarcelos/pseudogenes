process PREPARE_GENOME_FASTA {
    label 'light'
    publishDir "${params.outdir}/intermediate", mode: 'copy', pattern: 'genome.fna'

    input:
    path genome_input

    output:
    path 'genome.fna', emit: genome_fasta

    script:
    """
    if [[ "${genome_input}" == *.gz ]]; then
      gunzip -c "${genome_input}" > genome.fna
    else
      cp "${genome_input}" genome.fna
    fi
    """
}
