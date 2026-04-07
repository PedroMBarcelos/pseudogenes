process DOWNLOAD_GENOME {
    label 'light'
    publishDir "${params.outdir}/references", mode: 'copy'

    input:
    val url

    output:
    path 'genome.fna.gz'

    script:
    """
    wget -O genome.fna.gz "${url}"
    """
}
