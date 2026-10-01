import json
import pandas as pd
from pathlib import Path
from typing import List

from models.models import BlastHit, GenomicRegion, PseudogeneAnnotation


BLAST_COLUMNS = [
    "qseqid",
    "qlen",
    "sseqid",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "qseq",
    "sseq",
    "evalue",
    "bitscore",
    "score",
    "length",
    "pident",
    "nident",
    "mismatch",
    "positive",
    "gapopen",
    "gaps",
    "ppos",
    "sframe",
    "sstrand",
    "qcovs",
    "qcovhsp",
]


def select_best_protein(proteins_in_region: dict):
    best_protein = None
    max_density = 0.0
    for qseqid, data in proteins_in_region.items():
        if data['total_length'] > 0:
            density = data['total_score'] / data['total_length']
            if density > max_density:
                max_density = density
                best_protein = qseqid

    return best_protein


def parse_merged_bed(merged_bed_path: Path) -> List[GenomicRegion]:
    """
    Parse merged genomic regions TSV file.

    Expected columns:
        RegionID
        Chromosome
        Start
        End
        Strand
        ParentProtein
        ScoreDensity
        NumHSPs
        TotalFrameshifts
        InternalStops
        TotalDisablements

    Args:
        merged_bed_path: Path to merged genomic regions TSV file.

    Returns:
        List of GenomicRegion objects.
    """

    if not merged_bed_path.exists():
        raise FileNotFoundError(
            f"Arquivo de regiões genômicas não encontrado: {merged_bed_path}"
        )

    df = pd.read_csv(
        merged_bed_path,
        sep="\t",
        header=0,
    )

    regions = []

    for _, row in df.iterrows():
        region = GenomicRegion(
            region_id=str(row["RegionID"]),
            chromosome=str(row["Chromosome"]),
            start=int(row["Start"]),
            end=int(row["End"]),
            strand=str(row["Strand"]),
            parent_protein=str(row["ParentProtein"]),
            score_density=float(row["ScoreDensity"]),
            num_hsps=int(row["NumHSPs"]),
            total_frameshifts=int(row["TotalFrameshifts"]),
            internal_stops=int(row["InternalStops"]),
            total_disablements=int(row["TotalDisablements"]),
        )

        regions.append(region)

    return regions


def load_blast_hits(blast_path: Path) -> List[BlastHit]:
    """
    Load BLAST tblastn output and convert it to BlastHit objects.

    The BLAST output is expected to be tabular format (outfmt 6)
    without a header.
    """

    if not blast_path.exists():
        raise FileNotFoundError(
            f"Arquivo de output do BLAST não encontrado: {blast_path}"
        )

    df = pd.read_csv(
        blast_path,
        sep="\t",
        header=None,
        names=BLAST_COLUMNS,
    )

    hits = [
        BlastHit(**row)
        for row in df.to_dict("records")
    ]

    return hits


def load_pseudogene_annotations_from_json(json_path: Path) -> List[PseudogeneAnnotation]:
    """Carrega PseudogeneAnnotations de um arquivo JSON Lines."""
    annotations = []
    with json_path.open('r') as f:
        for line in f:
            line = line.strip()
            if line:
                d = json.loads(line)
                annotations.append(PseudogeneAnnotation.from_dict(d))
    return annotations