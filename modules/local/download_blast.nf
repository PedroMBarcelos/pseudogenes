process DOWNLOAD_BLAST {
    label 'light'
    publishDir "${projectDir}/", mode: 'copy'

    input:
    val url

    output:
    path 'ncbi-blast-2.17.0+/bin'

    script:
    """
    mkdir -p ncbi-blast-2.17.0+
    wget -O ncbi-blast-2.17.0+-x64-linux.tar.gz "${url}"
    
    tar -xvzf ncbi-blast-2.17.0+-x64-linux.tar.gz -C ncbi-blast-2.17.0+ --strip-components=1
    
    rm -rf fasta-36.3.8i-linux64.tar.gz
    """
}
