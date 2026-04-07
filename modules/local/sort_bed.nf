process SORT_BED {
    label 'light'
    publishDir "${params.outdir}/bed", mode: 'copy', pattern: 'protein_hits.sorted.bed'

    input:
    path bed_hits

    output:
    path 'protein_hits.sorted.bed', emit: sorted_bed

    script:
    """
    sort -k1,1 -k2,2n "${bed_hits}" > protein_hits.sorted.bed
    """
}
