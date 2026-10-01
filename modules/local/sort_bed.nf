process SORT_BED {
    tag { "${meta.id}" }
    label 'light'
    publishDir({ "${params.outdir}/${meta.id}/bed" }, mode: 'copy', pattern: '*.sorted.bed')

    input:
    tuple val(meta), path(bed_hits)

    output:
    tuple val(meta), path("${meta.id}.protein_hits.sorted.bed"), emit: sorted_bed

    script:
    """
    sort -k1,1 -k2,2n "${bed_hits}" > ${meta.id}.protein_hits.sorted.bed
    """
}