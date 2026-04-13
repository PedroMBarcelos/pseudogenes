process BLASTP_NULL_MODEL {
    label 'heavy'
    publishDir "${params.outdir}/blast", mode: 'copy', pattern: '*.blastp.out'

    input:
    path query_fasta
    val db_prefix
    path db_files

    output:
    path 'uniprot_sprot_shuffle.blastp.out', emit: blastp_out

    script:
    """
    "${params.blast_bin_dir}/blastp" \
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
      -evalue ${params.blast_evalue_null} \
      -max_target_seqs ${params.blast_max_target_seqs} \
      -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovs qcovhsp" \
      -query "${query_fasta}" \
      -db "${db_prefix}" \
      -out uniprot_sprot_shuffle.blastp.out
    """
}
