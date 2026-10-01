import logging
import pandas as pd
from pathlib import Path
from typing import List
from datetime import datetime

from models.models import PseudogeneAnnotation
from utils import load_pseudogene_annotations_from_json
from pseudogene_detection import calculate_subject_coverage_nt

logger = logging.getLogger(__name__)


def _build_gff3_attributes(
    annotation: PseudogeneAnnotation,
    index: int
) -> str:
    """
    Build attributes column (column 9) for GFF3.
    
    Format: key1=value1;key2=value2;...
    
    Args:
        annotation: Pseudogene annotation
        index: Annotation index (for unique ID)
    
    Returns:
        String with formatted attributes
    """
    attributes = []
    
    # Unique feature ID
    feature_id = f"lysozyme_{index}"
    attributes.append(f"ID={feature_id}")
    
    
    # Functional status
    status = "Pseudogene" if annotation.is_pseudogene else "Functional"
    attributes.append(f"Status={status}")
    
    # Small ORF status
    if annotation.is_small_orf:
        attributes.append("Small_ORF=True")
    
    # Reference protein
    best_protein_id = annotation.region_annotation.best_protein.protein_id
    if best_protein_id:
        attributes.append(f"Ref_Protein={best_protein_id}")
    
    # Detected mutation types
    mutations = []
    if annotation.disablements.non_synonymous_substitutions > 0:
        mutations.append(f"NonSynSubst:{annotation.disablements.non_synonymous_substitutions}")
    if annotation.disablements.in_frame_indels > 0:
        mutations.append(f"InFrameIndels:{annotation.disablements.in_frame_indels}")
    if annotation.disablements.frameshifts > 0:
        mutations.append(f"Frameshifts:{annotation.disablements.frameshifts}")
    if annotation.disablements.premature_stop_codons > 0:
        mutations.append(f"PrematureStops:{annotation.disablements.premature_stop_codons}")
    if annotation.disablements.missing_start_codon > 0:
        mutations.append(f"MissingStart:{annotation.disablements.missing_start_codon}")
    if annotation.disablements.missing_stop_codon > 0:
        mutations.append(f"MissingStop:{annotation.disablements.missing_stop_codon}")
    
    if mutations:
        attributes.append(f"Mutations={','.join(mutations)}")
    
    # Total inactivating mutations
    attributes.append(f"Total_Disablements={annotation.disablements.total_disablements}")
    
    # Score density
    score_density = annotation.region_annotation.best_protein.score_density
    attributes.append(f"Score_Density={score_density:.2f}")
    
    # Descriptive name
    name = f"Lysozyme_homolog_{index}"
    if annotation.is_pseudogene:
        name += "_pseudogene"
    attributes.append(f"Name={name}")
    
    return ";".join(attributes)


def save_pseudogene_annotations(
    annotations: List[PseudogeneAnnotation],
    output_path
) -> None:
    """
    Salva anotações de pseudogenes em arquivo TSV.
    
    Args:
        annotations: Lista de anotações de pseudogenes
        output_path: Caminho para o arquivo de saída
    """
    
    data = [ann.to_tsv_dict() for ann in annotations]
    df = pd.DataFrame(data)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)


def export_to_gff3(
    pseudogene_annotations: List[PseudogeneAnnotation],
    output_path: Path,
    source: str = "pseudogenes_pipeline"
) -> None:
    """
    Export pseudogene annotations to GFF3 format.
    
    GFF3 format:
    seqid source type start end score strand phase attributes
    
    Args:
        pseudogene_annotations: List of pseudogene annotations
        genome_id: Genome identifier
        output_path: Path to output GFF3 file
        source: Tool/pipeline name (column 2)
    """
    logger.debug(f"Exporting {len(pseudogene_annotations)} annotations to GFF3: {output_path}")
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        # Write GFF3 header
        f.write("##gff-version 3\n")
        f.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        f.write(f"##source lysozyme_pipeline v1.0\n")
        
        # Write each annotation
        for i, annotation in enumerate(pseudogene_annotations, 1):
            # Determine feature type
            feature_type = "pseudogene" if annotation.is_pseudogene else "CDS"
            
            # Determine strand (assume + if no information)
            strand = "+"
            
            # Determine score (use best protein score_density)
            score_density = annotation.region_annotation.best_protein.score_density
            score = f"{score_density:.2f}" if score_density > 0 else "."
            
            # Build attributes column (9)
            attributes = _build_gff3_attributes(annotation, i)
            
            # GFF3 line
            gff3_line = "\t".join([
                annotation.region_annotation.region.chromosome,  # seqid
                source,                                   # source
                feature_type,                             # type
                str(annotation.region_annotation.region.start),  # start
                str(annotation.region_annotation.region.end),    # end
                score,                                    # score
                strand,                                   # strand
                ".",                                      # phase (not applicable)
                attributes                                # attributes
            ])
            
            f.write(gff3_line + "\n")
    
    logger.debug(f"GFF3 file saved: {output_path}")


def save_coverage_statistics(
    annotations: List[PseudogeneAnnotation],
    output_path
) -> None:
    """
    Save detailed coverage statistics for further analysis.
    
    Creates TSV with coverage ratio, alignment details, and classification
    for discussion about size-based validation approaches.
    
    Args:
        annotations: List of pseudogene annotations
        output_path: Path to coverage statistics output file
    """
    import pandas as pd
    
    logger.debug(f"Saving coverage statistics to: {output_path}")

    coverage_data = []
    for ann in annotations:
        hsps = ann.region_annotation.best_protein.hsps
        if hsps:
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            coverage_len = max_qend - min_qstart + 1
            ref_len = hsps[0].qlen
            coverage_ratio = coverage_len / ref_len if ref_len > 0 else 0
            
            # Calculate genomic length
            genomic_len_nt = ann.region_annotation.region.length
            genomic_len_aa = genomic_len_nt / 3
            
            # Calculate subject coverage (nt)
            alignment_coverage_nt = calculate_subject_coverage_nt(hsps)
            
            # Calculate reference coverage in nt (alignment span on reference * 3)
            reference_coverage_nt = coverage_len * 3
            
            # Calculate alignment/genomic ratio (User requested metric)
            # (reference_coverage_nt) / genomic_region_nt
            alignment_genomic_ratio = reference_coverage_nt / genomic_len_nt if genomic_len_nt > 0 else 0
            
            # Determine classification
            if ann.is_pseudogene:
                if ann.is_small_orf:
                    classification = 'Pseudogene (Small ORF)'
                else:
                    classification = 'Pseudogene (Detected)'
            else:
                if ann.is_small_orf:
                    classification = 'Functional (Small ORF)'
                else:
                    classification = 'Functional (Possible Gene)'
            
            coverage_data.append({
                'chromosome': ann.region_annotation.region.chromosome,
                'start': ann.region_annotation.region.start,
                'end': ann.region_annotation.region.end,
                'strand': ann.region_annotation.region.strand,
                'protein_id': ann.region_annotation.best_protein.protein_id,
                'is_pseudogene': ann.is_pseudogene,
                'is_small_orf': ann.is_small_orf,
                'classification': classification,
                'alignment_coverage_aa': coverage_len,
                'alignment_coverage_nt': alignment_coverage_nt,
                'reference_coverage_nt': reference_coverage_nt,
                'alignment_genomic_ratio': round(alignment_genomic_ratio, 3),
                'reference_length_aa': ref_len,
                'coverage_ratio': round(coverage_ratio, 3),
                'genomic_region_nt': genomic_len_nt,
                'genomic_region_aa': round(genomic_len_aa, 1),
                'num_hsps': len(hsps),
                'qstart_min': min_qstart,
                'qend_max': max_qend,
                'missing_start_codon': ann.disablements.missing_start_codon,
                'missing_stop_codon': ann.disablements.missing_stop_codon,
                'premature_stops': ann.disablements.premature_stop_codons,
                'frameshifts': ann.disablements.frameshifts
            })
    
    df = pd.DataFrame(coverage_data)
    
    if df.empty or 'coverage_ratio' not in df.columns:
        # Write an empty file or headers so downstream/next steps don't crash
        df.to_csv(output_path, sep='\t', index=False)
        return
    
    # Sort by coverage ratio (ascending) to highlight problematic cases
    df = df.sort_values('coverage_ratio', ascending=True)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    
    logger.debug(f"Coverage statistics saved: {len(coverage_data)} regions")


def generate_summary_report(annotations: List[PseudogeneAnnotation], output_path: Path, min_coverage: float = 0.8) -> str:
    """
    Gera relatório resumido das anotações.
    
    Args:
        annotations: Lista de anotações
        min_coverage: Limite de cobertura para relatório (default: 0.8)
    
    Returns:
        String com o relatório formatado
    """
    total_regions = len(annotations)
    
    # Classification Breakdown
    functional_anns = [ann for ann in annotations if not ann.is_pseudogene]
    pseudogene_anns = [ann for ann in annotations if ann.is_pseudogene]
    
    num_functional = len(functional_anns)
    num_pseudogenes = len(pseudogene_anns)
    
    # Functional Sub-categories
    func_possible_genes = sum(1 for ann in functional_anns if not ann.is_small_orf)
    func_small_orfs = sum(1 for ann in functional_anns if ann.is_small_orf)
    
    # Pseudogene Sub-categories
    pseudo_detected = sum(1 for ann in pseudogene_anns if not ann.is_small_orf)
    pseudo_small_orfs = sum(1 for ann in pseudogene_anns if ann.is_small_orf)
    
    # Estatísticas de mutações
    total_substitutions = sum(ann.disablements.non_synonymous_substitutions for ann in annotations)
    total_indels = sum(ann.disablements.in_frame_indels for ann in annotations)
    total_frameshifts = sum(ann.disablements.frameshifts for ann in annotations)
    total_missing_start = sum(ann.disablements.missing_start_codon for ann in annotations)
    total_missing_stop = sum(ann.disablements.missing_stop_codon for ann in annotations)
    total_premature_stops = sum(ann.disablements.premature_stop_codons for ann in annotations)
    
    # NOVA SEÇÃO: Análise de cobertura e tamanho
    coverage_stats = []
    for ann in annotations:
            
        hsps = ann.region_annotation.best_protein.hsps
        if hsps:
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            coverage_len = max_qend - min_qstart + 1
            ref_len = hsps[0].qlen
            coverage_ratio = coverage_len / ref_len if ref_len > 0 else 0
            
            coverage_stats.append({
                'coverage_len': coverage_len,
                'ref_len': ref_len,
                'ratio': coverage_ratio,
                'is_pseudogene': ann.is_pseudogene,
                'is_small_orf': ann.is_small_orf
            })
    
    # Helper function for stats
    def calc_stats(ratios):
        if not ratios:
            return 0, 0, 0, 0
        avg = sum(ratios) / len(ratios)
        mn = min(ratios)
        mx = max(ratios)
        below_threshold = sum(1 for r in ratios if r < min_coverage)
        return avg, mn, mx, below_threshold
    
    # DEBUG
    # print(f"DEBUG: min_coverage={min_coverage}")

    # 1. Functional - Possible Genes
    func_possible_ratios = [s['ratio'] for s in coverage_stats if not s['is_pseudogene'] and not s['is_small_orf']]
    fp_avg, fp_min, fp_max, fp_below = calc_stats(func_possible_ratios)
    fp_count = len(func_possible_ratios)

    # 2. Functional - Small ORFs
    func_small_ratios = [s['ratio'] for s in coverage_stats if not s['is_pseudogene'] and s['is_small_orf']]
    fs_avg, fs_min, fs_max, fs_below = calc_stats(func_small_ratios)
    fs_count = len(func_small_ratios)

    # 3. Pseudogenes - Detected
    pseudo_detected_ratios = [s['ratio'] for s in coverage_stats if s['is_pseudogene'] and not s['is_small_orf']]
    pd_avg, pd_min, pd_max, pd_below = calc_stats(pseudo_detected_ratios)
    pd_count = len(pseudo_detected_ratios)

    # 4. Pseudogenes - Small ORFs (FIXED variable references to use 'ps_')
    pseudo_small_ratios = [s['ratio'] for s in coverage_stats if s['is_pseudogene'] and s['is_small_orf']]
    ps_avg, ps_min, ps_max, ps_below = calc_stats(pseudo_small_ratios)
    ps_count = len(pseudo_small_ratios)
    
    report = f"""
╭──────────────────────────────────────────────────────────────────╮
│                   PSEUDOGENE ANNOTATION REPORT                   │
╰──────────────────────────────────────────────────────────────────╯

GENERAL SUMMARY:
  Total regions analyzed:          {total_regions}
  
  Classification:
    Functional Genes:                {num_functional} ({100*num_functional/total_regions if total_regions else 0:.1f}%)
      - Possible Genes:              {func_possible_genes} ({100*func_possible_genes/num_functional if num_functional else 0:.1f}%)
      - Small ORFs:                  {func_small_orfs} ({100*func_small_orfs/num_functional if num_functional else 0:.1f}%)
      
    Pseudogenes:                     {num_pseudogenes} ({100*num_pseudogenes/total_regions if total_regions else 0:.1f}%)
      - Detected Pseudogenes:        {pseudo_detected} ({100*pseudo_detected/num_pseudogenes if num_pseudogenes else 0:.1f}%)
      - Small ORFs:                  {pseudo_small_orfs} ({100*pseudo_small_orfs/num_pseudogenes if num_pseudogenes else 0:.1f}%)

MUTATION STATISTICS:
  Non-synonymous substitutions:    {total_substitutions}
  In-frame indels:                 {total_indels}
  Frameshifts:                     {total_frameshifts}
  Missing start codon:             {total_missing_start}
  Missing stop codon:              {total_missing_stop}
  Premature stop codons:           {total_premature_stops}
  
  Total inactivating mutations:    {total_substitutions + total_indels + total_frameshifts + total_missing_start + total_missing_stop + total_premature_stops}

REFERENCE PROTEIN COVERAGE ANALYSIS:
  (Ratio = Alignment Coverage / Reference Size)
  
  1. Functional - Possible Genes ({fp_count} regions):
    Mean coverage ratio:             {fp_avg:.2f}
    Minimum ratio:                   {fp_min:.2f}
    Maximum ratio:                   {fp_max:.2f}
    Regions with coverage <{int(min_coverage*100)}%:      {fp_below} ({100*fp_below/fp_count if fp_count else 0:.1f}%)

  2. Functional - Small ORFs ({fs_count} regions):
    Mean coverage ratio:             {fs_avg:.2f}
    Minimum ratio:                   {fs_min:.2f}
    Maximum ratio:                   {fs_max:.2f}
    Regions with coverage <{int(min_coverage*100)}%:      {fs_below} ({100*fs_below/fs_count if fs_count else 0:.1f}%)
  
  3. Pseudogenes - Detected ({pd_count} regions):
    Mean coverage ratio:             {pd_avg:.2f}
    Minimum ratio:                   {pd_min:.2f}
    Maximum ratio:                   {pd_max:.2f}
    Regions with coverage <{int(min_coverage*100)}%:      {pd_below} ({100*pd_below/pd_count if pd_count else 0:.1f}%)

  4. Pseudogenes - Small ORFs ({ps_count} regions):
    Mean coverage ratio:             {ps_avg:.2f}
    Minimum ratio:                   {ps_min:.2f}
    Maximum ratio:                   {ps_max:.2f}
    Regions with coverage <{int(min_coverage*100)}%:      {ps_below} ({100*ps_below/ps_count if ps_count else 0:.1f}%)


"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(report)
    
    return report


def write_final_files(pseudogene_annotations: List[PseudogeneAnnotation], output_dir: Path) -> None:
    pseudogenes_output = output_dir / "pseudogene_annotations_final.tsv"
    logger.info(f"Saving {len(pseudogene_annotations)} pseudogene annotations to: {pseudogenes_output}")
    save_pseudogene_annotations(pseudogene_annotations, pseudogenes_output)
    logger.info("Pseudogene annotations saved successfully")
    
    
    # Export GFF3
    gff3_output = output_dir / "pseudogene_annotations.gff3"
    export_to_gff3(pseudogene_annotations, gff3_output)
    
    
    # Save coverage statistics for detailed analysis
    coverage_output = output_dir / "coverage_statistics.tsv"
    save_coverage_statistics(pseudogene_annotations, coverage_output)
    
    
    # Generate summary report
    report_file = output_dir / "summary_report.txt"
    summary = generate_summary_report(pseudogene_annotations, report_file)
    print(summary)
    
    num_pseudogenes = sum(1 for ann in pseudogene_annotations if ann.is_pseudogene)
    logger.info(f"Complete: {len(pseudogene_annotations)} regions, {num_pseudogenes} pseudogenes")



import sys
if __name__ == "__main__":
    logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
        )
    pseudogene_annotations = load_pseudogene_annotations_from_json(Path(sys.argv[1]))
    write_final_files(pseudogene_annotations, Path(sys.argv[2]))