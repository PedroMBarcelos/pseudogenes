process TBLASTN_TO_BED {
    label 'light'
    publishDir "${params.outdir}/bed", mode: 'copy', pattern: 'protein_hits.bed'

    input:
    path filtered_tblastn

    output:
    path 'protein_hits.bed', emit: bed_hits

    script:
    """
    awk 'BEGIN {OFS="\\t"} {
      if (\$7 < \$8) {
        start = \$7 - 1; end = \$8; strand = "+";
      } else {
        start = \$8 - 1; end = \$7; strand = "-";
      }
      print \$3, start, end, \$1, \$12, strand, NR;
    }' "${filtered_tblastn}" > protein_hits.bed
    """
}
