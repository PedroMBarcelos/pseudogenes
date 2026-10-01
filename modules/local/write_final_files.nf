process WRITE_FINAL_FILES {
    publishDir "${params.outdir}/final", mode: 'copy'
    
    input:
        path pseudogene_annotations
    
    output:
        path "pseudogene_annotations_final.tsv"
        path "pseudogene_annotations.gff3"
        path "coverage_statistics.tsv"
        path "summary_report.txt"

    script:
        """
        python3 ${workflow.projectDir}/write_final_files.py ${pseudogene_annotations} ./
        """
}