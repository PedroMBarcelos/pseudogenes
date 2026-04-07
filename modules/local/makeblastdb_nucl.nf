process MAKEBLASTDB_NUCL {
    label 'medium'
    publishDir "${params.outdir}/blastdb", mode: 'copy'

    input:
    path fasta

    output:
    val 'ecolik12db', emit: db_prefix
    path 'ecolik12db.*', emit: db_files

    script:
    """
    "${params.blast_bin_dir}/makeblastdb" -in "${fasta}" -dbtype nucl -out ecolik12db
    """
}
