include { DOWNLOAD_UNIPROT } from '../modules/local/download_uniprot'
include { DOWNLOAD_GENOME } from '../modules/local/download_genome'
include { DOWNLOAD_BLAST } from '../modules/local/download_blast'
include { DOWNLOAD_FASTA } from '../modules/local/download_fasta'
include { DOWNLOAD_GFF } from '../modules/local/download_gff'
include { PREPARE_UNIPROT_AND_SHUFFLE } from '../modules/local/prepare_uniprot_and_shuffle'
include { PREPARE_GENOME_FASTA } from '../modules/local/prepare_genome_fasta_multiple'
include { MAKEBLASTDB_PROT } from '../modules/local/makeblastdb_prot'
include { MAKEBLASTDB_NUCL } from '../modules/local/makeblastdb_nucl_multiple'
include { BLASTP_NULL_MODEL } from '../modules/local/blastp_null_model'
include { EXTRACT_MIN_EVALUE } from '../modules/local/extract_min_evalue'
include { TBLASTN_GENOME } from '../modules/local/tblastn_genome_multiple'
include { FILTER_TBLASTN } from '../modules/local/filter_tblastn_multiple'
include { EXTRACT_FRAGMENT_CHUNKS } from '../modules/local/extract_fragment_chunks_multiple'
include { SSEARCH_REALIGN_CHUNK } from '../modules/local/ssearch_realign_chunk_multiple'
include { CONCAT_SSEARCH_RESULTS } from '../modules/local/concat_ssearch_results_multiple'
include { MAKE_MIN_EVALUE_FILE } from '../modules/local/make_min_evalue_file'
include { TBLASTN_TO_BED } from '../modules/local/tblastn_to_bed_multiple'
include { SORT_BED } from '../modules/local/sort_bed_multiple'
include { MERGE_BED } from '../modules/local/merge_bed_multiple'
include { ANALYZE_REGIONS } from '../modules/local/analyze_regions_multiple'
include { ANNOTATE_REGIONS_WITH_BEST_PROTEINS } from "../modules/local/annotate_regions_proteins_multiple"
include { APPLY_FINAL_FILTERS } from "../modules/local/apply_final_filters_multiple"
include { PSEUDOGENE_DETECTION } from "../modules/local/pseudogene_detection_multiple"
include { WRITE_FINAL_FILES } from "../modules/local/write_final_files_multiple"

workflow PSEUDOGENES_TEST {
    main:
    
    // 1. Global Shared Resources (Downloaded/Prepared once)
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

    if (params.run_null_model) {
        PREPARE_UNIPROT_AND_SHUFFLE(uniprot_source)
        MAKEBLASTDB_PROT(PREPARE_UNIPROT_AND_SHUFFLE.out.shuffled_fasta, blast_bin_dir)
    }

    def min_evalue_ch
    if (params.min_evalue_file) {
        min_evalue_ch = Channel.value(file(params.min_evalue_file))
    } else if (params.min_evalue_manual != null) {
        MAKE_MIN_EVALUE_FILE(Channel.value(params.min_evalue_manual))
        min_evalue_ch = MAKE_MIN_EVALUE_FILE.out.min_evalue.first()
    } else if (params.run_null_model) {
        BLASTP_NULL_MODEL(PREPARE_UNIPROT_AND_SHUFFLE.out.uniprot_fasta, MAKEBLASTDB_PROT.out.db_prefix, MAKEBLASTDB_PROT.out.db_files)
        EXTRACT_MIN_EVALUE(BLASTP_NULL_MODEL.out.blastp_out)
        min_evalue_ch = EXTRACT_MIN_EVALUE.out.min_evalue.first()
    } else {
        error "No threshold source configured. Set --min_evalue_file, --min_evalue_manual, or enable --run_null_model."
    }
    
    // 2. Multi-Genome Input Channel with Safe Native CDS Counting
    def min_cds_threshold = 500 

    def genomes_ch = Channel
        .fromPath("${params.genomes_dir}/GCA_*", type: 'dir')
        .mix(Channel.fromPath("${params.genomes_dir}/GCF_*", type: 'dir'))
        .map { dir ->
            def id = dir.getName()
            def fastaList = files("${dir}/*_genomic.fna")
            def fasta = fastaList instanceof Collection ? fastaList[0] : (fastaList.exists() ? fastaList : null)
            def gffList = files("${dir}/*.gff*")
            def gffFile = gffList instanceof Collection ? gffList[0] : (gffList.exists() ? gffList : null)

            int cdsCount = 0
            if (gffFile && gffFile.exists()) {
                // Safely count CDS lines natively in Groovy without external command failures
                gffFile.eachLine { line ->
                    if (!line.startsWith("#")) {
                        def parts = line.tokenize('\t')
                        if (parts.size() > 2 && parts[2] == 'CDS') {
                            cdsCount += 1
                        }
                    }
                }
            }

            return [ [id: id, cds_count: cdsCount], fasta, gffFile ]
        }
        .filter { meta, fasta, gff -> 
            if (fasta == null || gff == null) {
                return false
            }
            if (meta.cds_count < min_cds_threshold) {
                log.info "Skipping ${meta.id}: Low annotation quality (${meta.cds_count} CDS features)."
                return false
            }
            // log.info "Keeping well-annotated genome ${meta.id} (${meta.cds_count} CDS features)."
            return true
        }

    // 3. Pass the complete bundle straight into preparation
    PREPARE_GENOME_FASTA(genomes_ch)

    // 4. Feed the output directly into MAKEBLASTDB_NUCL without any joins!
    MAKEBLASTDB_NUCL(
        PREPARE_GENOME_FASTA.out.genome_and_gff, 
        //blast_bin_dir
    )

    // 4. Downstream Pipeline Execution per Genome
    
    // Select the correct uniprot channel (shuffled or standard file wrapped as a channel)
    def uniprot_ch = params.run_null_model ? 
        PREPARE_UNIPROT_AND_SHUFFLE.out.uniprot_fasta : 
        Channel.value(files(params.uniprot_fasta))

    // Combine with the nucleotide genome database bundle
    def tblastn_input_ch = uniprot_ch
        .combine(MAKEBLASTDB_NUCL.out.db_bundle)
        .map { uniprot, meta, prefix, files, masked_fasta -> 
            [ meta, uniprot, prefix, files ] 
        }

    TBLASTN_GENOME(tblastn_input_ch)
    
    // 1. Filter tblastn output
    FILTER_TBLASTN(TBLASTN_GENOME.out.tblastn_out, min_evalue_ch)

    // 2. Extract chunks per genome
    EXTRACT_FRAGMENT_CHUNKS(FILTER_TBLASTN.out.filtered_tblastn)

    // 3. Transpose so every chunk file runs in parallel, retaining the meta tracking
    SSEARCH_REALIGN_CHUNK(
        EXTRACT_FRAGMENT_CHUNKS.out.fragment_chunk.transpose()
    )
    
    CONCAT_SSEARCH_RESULTS(
        SSEARCH_REALIGN_CHUNK.out.chunk_result.groupTuple()
    )

    TBLASTN_TO_BED(FILTER_TBLASTN.out.filtered_tblastn)
    SORT_BED(TBLASTN_TO_BED.out.bed_hits)
    MERGE_BED(SORT_BED.out.sorted_bed)

    // Join the three tuple streams together using the meta key
    def analysis_input_ch = CONCAT_SSEARCH_RESULTS.out.ssearch_results
        .join(SORT_BED.out.sorted_bed)
        .join(MERGE_BED.out.merged_bed)

    // Run region analysis per genome
    ANALYZE_REGIONS(analysis_input_ch, min_evalue_ch)

    // Join the analysis report and the SSEARCH/HSP results streams together
    def annotation_input_ch = ANALYZE_REGIONS.out.final_report
        .join(CONCAT_SSEARCH_RESULTS.out.ssearch_results)

    // Run region annotation per genome
    ANNOTATE_REGIONS_WITH_BEST_PROTEINS(
        ANALYZE_REGIONS.out.final_report
            .join(FILTER_TBLASTN.out.filtered_tblastn)
    )

    // Join the annotation output stream with your genome FASTA stream by meta.id
    def pseudogene_input_ch = ANNOTATE_REGIONS_WITH_BEST_PROTEINS.out.annotations
        .join(MAKEBLASTDB_NUCL.out.db_bundle.map { meta, db, db_files, masked_fasta -> [meta, masked_fasta] })

    // Run pseudogene detection per genome
    PSEUDOGENE_DETECTION(pseudogene_input_ch)

    // Run final filters per genome
    APPLY_FINAL_FILTERS(
        PSEUDOGENE_DETECTION.out.pseudogenes,
        params.min_coverage,
        params.final_min_identity
    )
    
    // Generate final output files per genome
    WRITE_FINAL_FILES(APPLY_FINAL_FILTERS.out.filtered_annotations)

    emit:
    report = ANALYZE_REGIONS.out.final_report
}