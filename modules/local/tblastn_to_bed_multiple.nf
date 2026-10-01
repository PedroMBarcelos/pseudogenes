process TBLASTN_TO_BED {
    tag { "${meta.id}" }
    label 'light'
    publishDir({ "${params.outdir}/${meta.id}/bed" }, mode: 'copy', pattern: '*.bed')

    input:
    tuple val(meta), path(filtered_tblastn)

    output:
    tuple val(meta), path("${meta.id}.protein_hits.bed"), emit: bed_hits

    script:
    """
    awk 'BEGIN {OFS="\\t"} {
      if (\$7 < \$8) {
        start = \$7 - 1; end = \$8; strand = "+";
      } else {
        start = \$8 - 1; end = \$7; strand = "-";
      }
      print \$3, start, end, \$1, \$12, strand, NR;
    }' "${filtered_tblastn}" > ${meta.id}.protein_hits.bed
    """
}