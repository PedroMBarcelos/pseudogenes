process APPLY_FINAL_FILTERS {
    publishDir "${params.outdir}/final", mode: 'copy'
    
    input:
        path pseudogene_annotations
        val min_coverage
        val final_min_identity

    output:
        path "pseudogene_annotations_with_coverage.jsonl"

    script:
        """
        python3 ${workflow.projectDir}/apply_final_filters.py ${pseudogene_annotations} ${min_coverage} ${final_min_identity} ./
        """ 

}