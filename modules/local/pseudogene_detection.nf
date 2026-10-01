process PSEUDOGENE_DETECTION {
    input:
        path region_annotations
        path genome_fasta
    
    output:
        path "initial_pseudogene_annotations.jsonl"

    script:
        """
        python3 ${workflow.projectDir}/pseudogene_detection.py ${region_annotations} ${genome_fasta} initial_pseudogene_annotations.jsonl
        """
}