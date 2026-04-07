process DOWNLOAD_UNIPROT {
    label 'light'
    publishDir "${params.outdir}/references", mode: 'copy'

    input:
    val url

    output:
    path 'uniprot_sprot.fasta.gz'

    script:
    """
    wget -O uniprot_sprot.fasta.gz "${url}"
    """
}
