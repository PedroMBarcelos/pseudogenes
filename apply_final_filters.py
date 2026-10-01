import json

import logging
from pathlib import Path
import sys
from typing import List

from models.models import PseudogeneAnnotation

logger = logging.getLogger(__name__)

def apply_coverage(pseudogene_annotations: List[PseudogeneAnnotation], min_coverage: float, output_path: Path) -> List[PseudogeneAnnotation]:
    filtered_annotations = []
        
    for ann in pseudogene_annotations:
        hsps = ann.region_annotation.best_protein.hsps
        if hsps:
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            coverage_len = max_qend - min_qstart + 1
            ref_len = hsps[0].qlen
            coverage_ratio = coverage_len / ref_len if ref_len > 0 else 0
            
            if coverage_ratio >= min_coverage:
                filtered_annotations.append(ann)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        for ann in filtered_annotations:
            f.write(json.dumps(ann.to_json_dict()) + '\n')  # JSON Lines: 1 objeto por linha
    
    return filtered_annotations


def apply_final_identity(pseudogene_annotations: List[PseudogeneAnnotation], final_min_identity: float,
                         output_path: Path) -> List[PseudogeneAnnotation]:
    threshold_pct = final_min_identity * 100 if final_min_identity <= 1.0 else final_min_identity
    
    filtered_annotations = []
    for ann in pseudogene_annotations:
        if max(hsp.pident for hsp in ann.region_annotation.best_protein.hsps) >= threshold_pct:
            filtered_annotations.append(ann)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
         for ann in filtered_annotations:
            f.write(json.dumps(ann.to_json_dict()) + '\n')  # JSON Lines: 1 objeto por linha

    return filtered_annotations


def apply_final_filters(pseudogene_annotations: List[PseudogeneAnnotation], min_coverage: float,
                        final_min_identity: float, output_dir: Path):
    """Apply coverage and final identity filters."""
    # --- Coverage Filter ---
    if min_coverage > 0:
        logger.info(
            f"Applying coverage filter: >= {min_coverage * 100:.1f}%"
        )

        original_count = len(pseudogene_annotations)

        annotations_with_coverage_path = (
            output_dir / "pseudogene_annotations_with_coverage.jsonl"
        )

        pseudogene_annotations = apply_coverage(
            pseudogene_annotations,
            min_coverage,
            annotations_with_coverage_path,
        )

        logger.info(
            f"  Filtered "
            f"{original_count - len(pseudogene_annotations)} regions. "
            f"Remaining: {len(pseudogene_annotations)}"
        )

    # --- Final Identity Filter ---
    if final_min_identity > 0:
        threshold_pct = (
            final_min_identity * 100
            if final_min_identity <= 1.0
            else final_min_identity
        )

        logger.info(
            f"Applying final identity filter: >= {threshold_pct:.1f}%"
        )

        original_count = len(pseudogene_annotations)

        annotations_with_final_identity_path = (
            output_dir / "pseudogene_annotations_with_final_identity.jsonl"
        )

        pseudogene_annotations = apply_final_identity(
            pseudogene_annotations,
            final_min_identity,
            annotations_with_final_identity_path,
        )

        logger.info(
            f"  Filtered "
            f"{original_count - len(pseudogene_annotations)} regions. "
            f"Remaining: {len(pseudogene_annotations)}"
        )

    return pseudogene_annotations

from utils import load_pseudogene_annotations_from_json

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    pseudogene_annotations = load_pseudogene_annotations_from_json(Path(sys.argv[1]))
    apply_final_filters(pseudogene_annotations, float(sys.argv[2]), float(sys.argv[3]), Path(sys.argv[4]))