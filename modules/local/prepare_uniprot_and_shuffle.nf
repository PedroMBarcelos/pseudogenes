process PREPARE_UNIPROT_AND_SHUFFLE {
    label 'medium'
    publishDir "${params.outdir}/intermediate", mode: 'copy', pattern: 'uniprot_sprot.fasta*'

    input:
    path uniprot_input

    output:
    path 'uniprot_sprot.fasta', emit: uniprot_fasta
    path 'uniprot_sprot.fasta.shuffle', emit: shuffled_fasta

    script:
    """
    if [[ "${uniprot_input}" == *.gz ]]; then
      cp "${uniprot_input}" uniprot_input.fasta.gz
    else
      gzip -c "${uniprot_input}" > uniprot_input.fasta.gz
    fi

    python3 "${projectDir}/Python_translation.py" uniprot_input.fasta.gz

    gunzip -c uniprot_input.fasta.gz > uniprot_sprot.fasta
    cp uniprot_input.fasta.gz.shuffle uniprot_sprot.fasta.shuffle
    """
}
