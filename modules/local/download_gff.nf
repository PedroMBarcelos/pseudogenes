process DOWNLOAD_GFF {
    label 'light'
    publishDir "${params.outdir}/references", mode: 'copy'

    input:
    val url

    output:
    path 'genome.gff'

    script:
    """
    wget -O genome.gff.gz "${url}"
    gunzip genome.gff.gz
    """
}
