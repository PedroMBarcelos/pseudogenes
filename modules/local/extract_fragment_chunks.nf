process EXTRACT_FRAGMENT_CHUNKS {
    label 'light'
    publishDir "${params.outdir}/intermediate", mode: 'copy', pattern: 'frag_chunk_*.tsv'

    input:
    path filtered_tblastn

    output:
    path 'frag_chunk_*.tsv', emit: fragment_chunk

    script:
    """
    awk 'BEGIN{OFS="\t"} {print NR, \$9, \$10}' "${filtered_tblastn}" > fragments_with_index.tsv
    
    # Ensure at least one output file is created, even if input is empty
    if [ -s fragments_with_index.tsv ]; then
        split -d -l ${params.ssearch_chunk_size} --additional-suffix=.tsv fragments_with_index.tsv frag_chunk_
    else
        touch frag_chunk_00.tsv
    fi
    """
}
