process MAKEBLASTDB_NUCL {
    label 'medium'
    publishDir({ "${params.outdir}/blastdb" }, mode: 'copy')

    input:
    path fasta
    path blast_bin_dir
    path gff_source

    output:
    val 'ecolik12db', emit: db_prefix
    path 'ecolik12db.*', emit: db_files

    script:
    """
    sed -i 's/\r//g' "${gff_source}"
    awk '$3=="CDS" {print $1 "\t" $4-1 "\t" $5}' "${gff_source}" | sort -k1,1 -k2,2n > cds_regions.bed
    bedtools maskfasta -fi "${projectDir}/results_nf/intermediate/genome.fna" -bed cds_regions.bed -fo "${projectDir}/results_nf/intermediate/genoma_cds_mascarado.fasta"
    "${blast_bin_dir}/makeblastdb" -in "${projectDir}/results_nf/intermediate/genoma_cds_mascarado.fasta" -dbtype nucl -out ecolik12db

    """
}
