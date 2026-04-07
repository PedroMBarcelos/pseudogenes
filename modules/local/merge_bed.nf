process MERGE_BED {
    label 'light'
    publishDir "${params.outdir}/bed", mode: 'copy', pattern: 'merged_protein_features.bed'

    input:
    path sorted_bed

    output:
    path 'merged_protein_features.bed', emit: merged_bed

    script:
    """
    bedtools merge -s -i "${sorted_bed}" -c 6 -o first > merged_protein_features.bed
    """
}
