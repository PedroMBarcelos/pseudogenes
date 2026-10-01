process ANNOTATE_REGIONS_WITH_BEST_PROTEINS {
    publishDir({ "${params.outdir}" }, mode: 'copy')

    input:
        path regions_path
        path hsps_path
    
    output:
    path "region_annotations.jsonl"

    script:
    """
    python3 ${workflow.projectDir}/annotate_regions_with_best_proteins.py ${regions_path} ${hsps_path} region_annotations.jsonl
    """
}