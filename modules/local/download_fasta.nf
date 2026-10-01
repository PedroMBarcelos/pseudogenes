process DOWNLOAD_FASTA {
    label 'light'
    publishDir({ "${projectDir}/" }, mode: 'copy')

    input:
    val url

    output:
    path 'fasta-36.3.8i/bin'

    script:
    """
    mkdir -p fasta-36.3.8i
    wget -O fasta-36.3.8i-linux64.tar.gz "${url}"
    
    tar -xvzf fasta-36.3.8i-linux64.tar.gz -C fasta-36.3.8i --strip-components=1
    
    rm -rf fasta-36.3.8i-linux64.tar.gz
    """
}
