process EXTRACT_FRAGMENT_CHUNKS {
    tag { "${meta.id}" }
    label 'light'
    publishDir({ "${params.outdir}/${meta.id}/intermediate" }, mode: 'copy', pattern: '${meta.id}_frag_chunk_*.tsv')

    input:
    tuple val(meta), path(filtered_tblastn)

    output:
    tuple val(meta), path("${meta.id}_frag_chunk_*.tsv"), emit: fragment_chunk

    script:
    """
    awk 'BEGIN{OFS="\t"} {print NR, \$9, \$10}' "${filtered_tblastn}" > ${meta.id}_fragments_with_index.tsv
    
    # Ensure at least one output file is created safely with unique prefixes
    if [ -s ${meta.id}_fragments_with_index.tsv ]; then
        split -d -l ${params.ssearch_chunk_size} --additional-suffix=.tsv ${meta.id}_fragments_with_index.tsv ${meta.id}_frag_chunk_
    else
        touch ${meta.id}_frag_chunk_00.tsv
    fi
    """
}