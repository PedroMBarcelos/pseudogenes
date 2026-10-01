include { DOWNLOAD_UNIPROT } from '../modules/local/download_uniprot'
include { DOWNLOAD_GENOME } from '../modules/local/download_genome'
include { DOWNLOAD_BLAST } from '../modules/local/download_blast'
include { DOWNLOAD_FASTA } from '../modules/local/download_fasta'
include { DOWNLOAD_GFF } from '../modules/local/download_gff'
include { PREPARE_UNIPROT_AND_SHUFFLE } from '../modules/local/prepare_uniprot_and_shuffle'
include { PREPARE_GENOME_FASTA } from '../modules/local/prepare_genome_fasta'
include { MAKEBLASTDB_PROT } from '../modules/local/makeblastdb_prot'
include { MAKEBLASTDB_NUCL } from '../modules/local/makeblastdb_nucl'
include { BLASTP_NULL_MODEL } from '../modules/local/blastp_null_model'
include { EXTRACT_MIN_EVALUE } from '../modules/local/extract_min_evalue'
include { TBLASTN_GENOME } from '../modules/local/tblastn_genome'
include { FILTER_TBLASTN } from '../modules/local/filter_tblastn'
include { EXTRACT_FRAGMENT_CHUNKS } from '../modules/local/extract_fragment_chunks'
include { SSEARCH_REALIGN_CHUNK } from '../modules/local/ssearch_realign_chunk'
include { CONCAT_SSEARCH_RESULTS } from '../modules/local/concat_ssearch_results'
include { MAKE_MIN_EVALUE_FILE } from '../modules/local/make_min_evalue_file'
include { TBLASTN_TO_BED } from '../modules/local/tblastn_to_bed'
include { SORT_BED } from '../modules/local/sort_bed'
include { MERGE_BED } from '../modules/local/merge_bed'
include { ANALYZE_REGIONS } from '../modules/local/analyze_regions'
include { ANNOTATE_REGIONS_WITH_BEST_PROTEINS } from "../modules/local/annotate_regions_proteins"
include { APPLY_FINAL_FILTERS } from "../modules/local/apply_final_filters"
include { PSEUDOGENE_DETECTION } from "../modules/local/pseudogene_detection"
include { WRITE_FINAL_FILES } from "../modules/local/write_final_files"

workflow PSEUDOGENES {
    main:
    def blast_bin_dir
    if (file(params.blast_bin_dir).exists()) {
        blast_bin_dir = Channel.of(file(params.blast_bin_dir))
    } else {
        blast_bin_dir = DOWNLOAD_BLAST(Channel.value(params.blast_url))
    }

    def fasta_bin_dir
    if (file(params.fasta_bin_dir).exists()) {
        fasta_bin_dir = Channel.of(file(params.fasta_bin_dir))
    } else {
        fasta_bin_dir = DOWNLOAD_FASTA(Channel.value(params.fasta_url))
    } 

    def uniprot_source
    if (params.uniprot_fasta && file(params.uniprot_fasta).exists()) {
        uniprot_source = Channel.of(file(params.uniprot_fasta))
    } else {
        uniprot_source = DOWNLOAD_UNIPROT(Channel.value(params.uniprot_url))
    }

    def genome_source
    if (params.genome_fasta && file(params.genome_fasta).exists()) {
        genome_source = Channel.of(file(params.genome_fasta))
    } else {
        genome_source = DOWNLOAD_GENOME(Channel.value(params.genome_url))
    }

    def gff_source
    if (params.genome_gff && file(params.genome_gff).exists()) {
        gff_source = Channel.of(file(params.genome_gff))
    } else {
        gff_source = DOWNLOAD_GFF(Channel.value(params.genome_gff_url))
    }

    PREPARE_UNIPROT_AND_SHUFFLE(uniprot_source)
    PREPARE_GENOME_FASTA(genome_source)


    MAKEBLASTDB_PROT(PREPARE_UNIPROT_AND_SHUFFLE.out.shuffled_fasta, blast_bin_dir)
    MAKEBLASTDB_NUCL(PREPARE_GENOME_FASTA.out.genome_fasta, blast_bin_dir, gff_source)

    def min_evalue_ch
    if (params.min_evalue_file) {
        min_evalue_ch = Channel.of(file(params.min_evalue_file))
    } else if (params.min_evalue_manual != null) {
        MAKE_MIN_EVALUE_FILE(Channel.value(params.min_evalue_manual))
        min_evalue_ch = MAKE_MIN_EVALUE_FILE.out.min_evalue
    } else if (params.run_null_model) {
        BLASTP_NULL_MODEL(PREPARE_UNIPROT_AND_SHUFFLE.out.uniprot_fasta, MAKEBLASTDB_PROT.out.db_prefix, MAKEBLASTDB_PROT.out.db_files)
        EXTRACT_MIN_EVALUE(BLASTP_NULL_MODEL.out.blastp_out)
        min_evalue_ch = EXTRACT_MIN_EVALUE.out.min_evalue
    } else {
        error "No threshold source configured. Set --min_evalue_file, --min_evalue_manual, or enable --run_null_model."
    }
    
    TBLASTN_GENOME(PREPARE_UNIPROT_AND_SHUFFLE.out.uniprot_fasta, MAKEBLASTDB_NUCL.out.db_prefix, MAKEBLASTDB_NUCL.out.db_files)
    FILTER_TBLASTN(TBLASTN_GENOME.out.tblastn_out, min_evalue_ch)

    EXTRACT_FRAGMENT_CHUNKS(FILTER_TBLASTN.out.filtered_tblastn)
    SSEARCH_REALIGN_CHUNK(EXTRACT_FRAGMENT_CHUNKS.out.fragment_chunk.flatten())
    CONCAT_SSEARCH_RESULTS(SSEARCH_REALIGN_CHUNK.out.chunk_result.collect())

    TBLASTN_TO_BED(FILTER_TBLASTN.out.filtered_tblastn)
    SORT_BED(TBLASTN_TO_BED.out.bed_hits)
    MERGE_BED(SORT_BED.out.sorted_bed)

    ANALYZE_REGIONS(
        CONCAT_SSEARCH_RESULTS.out.ssearch_results,
        SORT_BED.out.sorted_bed,
        MERGE_BED.out.merged_bed,
        min_evalue_ch
    )

    ANNOTATE_REGIONS_WITH_BEST_PROTEINS(
        ANALYZE_REGIONS.out.final_report,
        TBLASTN_GENOME.out.tblastn_out
    )

    PSEUDOGENE_DETECTION(
        ANNOTATE_REGIONS_WITH_BEST_PROTEINS.out,
        PREPARE_GENOME_FASTA.out.genome_fasta
    )

    APPLY_FINAL_FILTERS(
        PSEUDOGENE_DETECTION.out,
        params.min_coverage,
        params.final_min_identity
    )
    
    WRITE_FINAL_FILES(
        APPLY_FINAL_FILTERS.out
    )

    emit:
    report = ANALYZE_REGIONS.out.final_report
}
