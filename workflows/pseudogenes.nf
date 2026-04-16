include { DOWNLOAD_UNIPROT } from '../modules/local/download_uniprot'
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

workflow PSEUDOGENES {
    main:
    def uniprot_source
    if (params.uniprot_fasta && file(params.uniprot_fasta).exists()) {
        uniprot_source = Channel.of(file(params.uniprot_fasta))
    } else {
        uniprot_source = DOWNLOAD_UNIPROT(Channel.value(params.uniprot_url))
    }

    def genome_source
    if (!params.genome_fasta) {
        error "Missing required genome input. Provide --genome_fasta with a local FASTA (.fna/.fa) or gzipped FASTA (.gz) file."
    }
    if (!file(params.genome_fasta).exists()) {
        error "Genome file not found: ${params.genome_fasta}. Provide a valid local path via --genome_fasta."
    }
    genome_source = Channel.of(file(params.genome_fasta))

    PREPARE_UNIPROT_AND_SHUFFLE(uniprot_source)
    PREPARE_GENOME_FASTA(genome_source)

    MAKEBLASTDB_PROT(PREPARE_UNIPROT_AND_SHUFFLE.out.shuffled_fasta)
    MAKEBLASTDB_NUCL(PREPARE_GENOME_FASTA.out.genome_fasta)

    def min_evalue_ch
    if (params.run_null_model) {
        BLASTP_NULL_MODEL(PREPARE_UNIPROT_AND_SHUFFLE.out.uniprot_fasta, MAKEBLASTDB_PROT.out.db_prefix, MAKEBLASTDB_PROT.out.db_files)
        EXTRACT_MIN_EVALUE(BLASTP_NULL_MODEL.out.blastp_out)
        min_evalue_ch = EXTRACT_MIN_EVALUE.out.min_evalue
    } else {
        MAKE_MIN_EVALUE_FILE(Channel.value(params.min_evalue))
        min_evalue_ch = MAKE_MIN_EVALUE_FILE.out.min_evalue
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

    emit:
    report = ANALYZE_REGIONS.out.final_report
}
