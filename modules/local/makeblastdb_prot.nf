process MAKEBLASTDB_PROT {
    cache "lenient"
    label 'medium'
    publishDir({ "${params.outdir}/blastdb" }, mode: 'copy')

    input:
    path fasta
    path blast_bin_dir

    output:
    val 'uniprot_sprot_shuffle', emit: db_prefix
    path 'uniprot_sprot_shuffle.*', emit: db_files

    script:
    """
    "${blast_bin_dir}/makeblastdb" -in "${fasta}" -dbtype prot -out uniprot_sprot_shuffle
    """
}
