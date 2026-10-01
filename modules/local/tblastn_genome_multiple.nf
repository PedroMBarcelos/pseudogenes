process TBLASTN_GENOME {
    cache "lenient"
    tag "$meta.id"
    label 'medium'
    publishDir({ "${params.outdir}/${meta.id}/blast" }, mode: 'copy', pattern: '*.tblastn.out')

    input:
    tuple val(meta), path(query_fasta), val(db_prefix), path(db_files)

    output:
    tuple val(meta), path("${meta.id}.tblastn.out"), emit: tblastn_out

    script:
    """
    "${params.blast_bin_dir}/tblastn" \
      -num_threads ${task.cpus} \
      -word_size ${params.blast_word_size} \
      -gapopen ${params.blast_gapopen} \
      -gapextend ${params.blast_gapextend} \
      -matrix ${params.blast_matrix} \
      -threshold ${params.blast_threshold} \
      -comp_based_stats ${params.blast_comp_stats} \
      -seg yes \
      -soft_masking true \
      -lcase_masking \
      -evalue ${params.blast_evalue_tblastn} \
      -max_target_seqs ${params.blast_max_target_seqs} \
      -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sframe sstrand qcovs qcovhsp" \
      -query "${query_fasta}" \
      -db "${db_prefix}" \
      -out "${meta.id}.tblastn.out"
    """
}