process MAKEBLASTDB_NUCL {
    tag "$meta.id"
    label 'medium'
    publishDir "${params.outdir}/${meta.id}/blastdb", mode: 'copy', pattern: "${meta.id}_db.*"
    publishDir "${params.outdir}/${meta.id}/blastdb", mode: 'copy', pattern: "${meta.id}_genoma_cds_mascarado.fasta"

    input:
    tuple val(meta), path(fasta), path(gff_source)

    output:
    tuple val(meta), val("${meta.id}_db"), path("${meta.id}_db.*"), path("${meta.id}_genoma_cds_mascarado.fasta"), emit: db_bundle

    script:
    def db_name = "${meta.id}_db"
    """
    # 1. Extract CDS regions safely (ignoring comments and handling standard GFF3)
    awk '\$1 !~ /^#/ && \$3 == "CDS" {print \$1 "\t" \$4-1 "\t" \$5}' "${gff_source}" > ${meta.id}_cds_regions.bed
    
    # 2. Check if BED file actually got features; fallback gracefully if empty
    if [ ! -s ${meta.id}_cds_regions.bed ]; then
        echo "Warning: No CDS features found in ${gff_source} for ${meta.id}!"
    fi

    # 3. Mask fasta
    bedtools maskfasta -fi "${fasta}" -bed ${meta.id}_cds_regions.bed -fo ${meta.id}_genoma_cds_mascarado.fasta

    # 4. Build BLAST database
    "${params.blast_bin_dir}/makeblastdb" -in ${meta.id}_genoma_cds_mascarado.fasta -dbtype nucl -out ${db_name}
    """
}