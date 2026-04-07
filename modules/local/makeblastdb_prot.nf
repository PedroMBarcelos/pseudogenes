process MAKEBLASTDB_PROT {
    label 'medium'
    publishDir "${params.outdir}/blastdb", mode: 'copy'

    input:
    path fasta

    output:
    val 'uniprot_sprot_shuffle', emit: db_prefix
    path 'uniprot_sprot_shuffle.*', emit: db_files

    script:
    """
    "${params.blast_bin_dir}/makeblastdb" -in "${fasta}" -dbtype prot -out uniprot_sprot_shuffle
    """
}
