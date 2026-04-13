process SSEARCH_REALIGN_CHUNK {
    label 'heavy'
    publishDir "${params.outdir}/ssearch", mode: 'copy', pattern: '*.ssearch.out'

    input:
    path fragment_chunk

    output:
    path "${fragment_chunk.simpleName}.ssearch.out", emit: chunk_result

    script:
    """
    out_file="${fragment_chunk.simpleName}.ssearch.out"
    : > "\${out_file}"

    if command -v ssearch36 &> /dev/null; then
      ssearch_cmd="\$(command -v ssearch36)"
    elif [ -x "${params.ssearch_bin}" ]; then
      ssearch_cmd="${params.ssearch_bin}"
    else
      echo "WARNING: ssearch36 not found at ${params.ssearch_bin}" >&2
      echo "To install: download from https://fasta.bioch.virginia.edu/" >&2
      echo "For conda: conda install -c bioconda fasta" >&2
      exit 1
    fi

    while IFS=\$'\t' read -r line_num qseq_fragment sseq_fragment; do
      clean_qseq="\${qseq_fragment//-/}"
      clean_sseq="\${sseq_fragment//-/}"

      alignment_output=\$("\${ssearch_cmd}" \
        -b 1 -d 1 -k ${params.ssearch_k} -f ${params.ssearch_f} -g ${params.ssearch_g} \
        -m 10 -q -s ${params.ssearch_matrix} -z ${params.ssearch_z} -Z ${params.ssearch_Z} \
        <(echo -e ">query_frag_\${line_num}\\n\${clean_qseq}") \
        <(echo -e ">subject_frag_\${line_num}\\n\${clean_sseq}")
      )

      if grep -q "The best scores are:" <<< "\${alignment_output}"; then
        {
          echo "# Original TBLASTN hit line: \${line_num}"
          echo "\${alignment_output}"
        } >> "\${out_file}"
      fi
    done < "${fragment_chunk}"
    """
}
