import sys
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from utils import parse_merged_bed, load_blast_hits
from models.models import BlastHit, GenomicRegion, RegionAnnotation, ProteinHitGroup


def assign_hsps_to_regions(
    regions: List[GenomicRegion],
    all_hsps: List[BlastHit]
) -> Dict[Tuple[str, int, int], List[BlastHit]]:
    """
    Assign HSPs to merged genomic regions.
    
    Args:
        regions: List of merged genomic regions
        all_hsps: List of all HSPs
    
    Returns:
        Dictionary mapping (chromosome, start, end) to list of HSPs
    """
    
    region_hsps_map = defaultdict(list)
    
    for hsp in all_hsps:
        # HSP coordinates
        hsp_chrom = hsp.sseqid
        hsp_start = min(hsp.sstart, hsp.send)
        hsp_end = max(hsp.sstart, hsp.send)
        
        # Find regions that overlap with this HSP
        for region in regions:
            if region.chromosome != hsp_chrom:
                continue
            
            # Check overlap
            if hsp_start <= region.end and hsp_end >= region.start:
                region_key = (region.chromosome, region.start, region.end)
                region_hsps_map[region_key].append(hsp)
    
    return dict(region_hsps_map)


def group_hsps_by_protein(hsps: List[BlastHit]) -> Dict[str, List[BlastHit]]:
    """
    Group HSPs by query protein.
    
    Args:
        hsps: List of HSPs
    
    Returns:
        Dictionary mapping protein_id to list of its HSPs
    """
    protein_groups = defaultdict(list)

    for hsp in hsps:
        protein_groups[hsp.qseqid].append(hsp)
    
    return dict(protein_groups)


def select_best_protein_for_region(
    region: GenomicRegion,
    region_hsps: List[BlastHit]
) -> RegionAnnotation:
    """
    Select protein with highest score density for a region.
    
    Args:
        region: Merged genomic region
        region_hsps: List of HSPs that map to this region
    
    Returns:
        Region annotation with best protein selected
    """
    
    # Group HSPs by protein
    protein_groups = group_hsps_by_protein(region_hsps)

    # Calculate density for each protein
    protein_hit_groups = []
    for protein_id, hsps in protein_groups.items():
        hit_group = ProteinHitGroup.from_hsps(protein_id, hsps)
        protein_hit_groups.append(hit_group)


    # Select protein with highest density ("King of the Hill" algorithm)
    # Primary criterion: Highest Score Density
    # Tie-breaker: Highest Total Score (if densities are equal)
    best_protein = max(
        protein_hit_groups,
        key=lambda x: (x.score_density, x.total_score)
    )

     # Update query_ids list in region
    region.query_ids = [best_protein.protein_id]
    
    return RegionAnnotation(
        region=region,
        best_protein=best_protein,
        all_proteins=protein_hit_groups
    )


def annotate_regions_with_best_proteins(
    regions: List[GenomicRegion],
    all_hsps: List[BlastHit],
    output_path: Path,
):
    # Assign HSPs to regions
    region_hsps_map = assign_hsps_to_regions(regions, all_hsps)

    # For each region, select best protein
    annotations = []
    for region in regions:
        region_key = (region.chromosome, region.start, region.end)
        region_hsps = region_hsps_map.get(region_key, [])
        
        if not region_hsps:
            continue
        
        annotation = select_best_protein_for_region(region, region_hsps)
        annotations.append(annotation)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        for ann in annotations:
            f.write(json.dumps(ann.to_json_dict()) + '\n')  # JSON Lines: 1 objeto por linha

    return annotations
    

if __name__ == "__main__":
    genomic_region = parse_merged_bed(Path(sys.argv[1]))
    blast_hits = load_blast_hits(Path(sys.argv[2]))

    annotate_regions_with_best_proteins(genomic_region, blast_hits, Path(sys.argv[3]))