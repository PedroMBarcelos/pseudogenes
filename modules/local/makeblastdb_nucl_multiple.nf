process MAKEBLASTDB_NUCL {
    tag { "${meta.id}" }
    label 'medium'
    publishDir({ "${params.outdir}/${meta.id}/blastdb" }, mode: 'copy', pattern: "*_db.*")
    publishDir({ "${params.outdir}/${meta.id}/blastdb" }, mode: 'copy', pattern: "*_genoma_cds_mascarado.fasta")

    input:
    tuple val(meta), path(fasta), path(gff_source)

    output:
    tuple val(meta), val("${meta.id}_db"), path("${meta.id}_db.*"), path("${meta.id}_genoma_cds_mascarado.fasta"), emit: db_bundle

    script:
    def maskMode = params.mask_mode ?: 'protein_gene_preserve_pseudogene'
    def db_name = "${meta.id}_db"
    """
    # Masking implementation v1
    python3 "${projectDir}/scripts/build_genome_mask.py" \
        --fasta "${fasta}" \
        --gff "${gff_source}" \
        --mode "${maskMode}" \
        --bed "${meta.id}_mask.bed" \
        --report "${meta.id}_mask_summary.json"

    if [ "${maskMode}" = "none" ] || [ ! -s "${meta.id}_mask.bed" ]; then
        cp "${fasta}" "${meta.id}_genoma_cds_mascarado.fasta"
    else
        bedtools maskfasta \
            -fi "${fasta}" \
            -bed "${meta.id}_mask.bed" \
            -fo "${meta.id}_genoma_cds_mascarado.fasta"
    fi

    # 4. Build BLAST database
    "${params.blast_bin_dir}/makeblastdb" -in ${meta.id}_genoma_cds_mascarado.fasta -dbtype nucl -out ${db_name}
    """
}