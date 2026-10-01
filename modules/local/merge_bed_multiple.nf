process MERGE_BED {
    tag { "${meta.id}" }
    label 'light'
    publishDir({ "${params.outdir}/${meta.id}/bed" }, mode: 'copy', pattern: '*.merged.bed')

    input:
    tuple val(meta), path(sorted_bed)

    output:
    tuple val(meta), path("${meta.id}.merged_protein_features.bed"), emit: merged_bed

    script:
    """
    bedtools merge -s -i "${sorted_bed}" -c 6 -o first > ${meta.id}.merged_protein_features.bed
    """
}