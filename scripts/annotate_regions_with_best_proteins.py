import sys
import re
import json
import math
import argparse

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from utils import parse_merged_bed, load_blast_hits
from models.models import (
    BlastHit,
    GenomicRegion,
    RegionAnnotation,
    ProteinHitGroup,
)


# SSEARCH output must have a source-line marker followed
# by tabular m8-style alignments.
SOURCE_RE = re.compile(
    r"#\s*Original TBLASTN hit line:\s*(\d+)"
)


@dataclass(frozen=True)
class SsearchHit:
    source_line: int
    bitscore: float
    alignment_length: int
    evalue: float


def normalize_strand(value):
    if value is None:
        return None

    value = str(value).strip().lower()

    return {
        "+": "+",
        "plus": "+",
        "-": "-",
        "minus": "-",
    }.get(value)


def get_hsp_strand(hsp):
    strand = normalize_strand(
        getattr(hsp, "sstrand", None)
    )

    if strand is not None:
        return strand

    frame = getattr(hsp, "sframe", None)

    if frame is not None:
        frame = int(frame)
        if frame:
            return "+" if frame > 0 else "-"

    return normalize_strand(
        getattr(hsp, "strand", None)
    )


def load_indexed_tblastn(path):
    """
    Load filtered TBLASTN hits, preserving their
    original 1-based line numbers.

    The line numbers must match those used when
    generating SSEARCH fragments.
    """
    lines = path.read_text().splitlines()

    if any(not line.strip() or line.startswith("#")
           for line in lines):
        raise ValueError(
            "Filtered TBLASTN must contain only hit rows. "
            "Blank/comment lines would affect source indices."
        )

    hits = load_blast_hits(path)

    if len(hits) != len(lines):
        raise ValueError(
            f"TBLASTN parser returned {len(hits)} hits "
            f"for {len(lines)} lines. "
            "Cannot safely recover source indices."
        )

    return dict(enumerate(hits, start=1))


def parse_ssearch_results(path):
    """
    Parse concatenated FASTA36 SSEARCH -m 10 output.

    Each block begins with:
        # Original TBLASTN hit line: N

    Alignment statistics include:
        ; sw_bits: 31.1
        ; sw_expect: 0.0055
        ; sw_overlap: 37

    Returns:
        Dict[int, SsearchHit]

    When multiple alignments occur for the same original
    TBLASTN hit, retain the highest-scoring alignment.
    """

    best_hits = {}

    current_source = None
    current_alignment = {}

    def save_alignment():
        if current_source is None:
            return

        required = (
            "bitscore",
            "evalue",
            "alignment_length",
        )

        if not all(
            field in current_alignment
            for field in required
        ):
            return

        hit = SsearchHit(
            source_line=current_source,
            bitscore=current_alignment["bitscore"],
            alignment_length=current_alignment[
                "alignment_length"
            ],
            evalue=current_alignment["evalue"],
        )

        old = best_hits.get(current_source)

        if old is None or (
            hit.bitscore,
            -hit.evalue
        ) > (
            old.bitscore,
            -old.evalue
        ):
            best_hits[current_source] = hit

    with Path(path).open() as handle:
        for raw_line in handle:
            line = raw_line.strip()

            # Beginning of a new original TBLASTN hit
            match = SOURCE_RE.search(line)

            if match:
                save_alignment()

                current_source = int(match.group(1))
                current_alignment = {}
                continue

            if current_source is None:
                continue

            # Beginning of a new alignment record
            if line.startswith(">>>"):
                save_alignment()
                current_alignment = {}
                continue

            # FASTA36 -m 10 fields
            if line.startswith("; sw_bits:"):
                current_alignment["bitscore"] = float(
                    line.split(":", 1)[1].strip()
                )

            elif line.startswith("; sw_expect:"):
                current_alignment["evalue"] = float(
                    line.split(":", 1)[1].strip()
                )

            elif line.startswith("; sw_overlap:"):
                current_alignment[
                    "alignment_length"
                ] = int(
                    line.split(":", 1)[1].strip()
                )

    # Save the final alignment
    save_alignment()

    return best_hits


def assign_hsps_to_regions(
    regions: List[GenomicRegion],
    indexed_hsps: Dict[int, BlastHit],
):
    """
    Assign TBLASTN source-line indices to genomic regions.

    Assumes:
      - TBLASTN: 1-based inclusive coordinates
      - GenomicRegion: 0-based half-open BED coordinates

    Strand mismatches are rejected.
    """
    region_map = defaultdict(list)

    for source_line, hsp in indexed_hsps.items():
        hsp_start = min(hsp.sstart, hsp.send) - 1
        hsp_end = max(hsp.sstart, hsp.send)

        hsp_strand = get_hsp_strand(hsp)

        for idx, region in enumerate(regions):
            if region.chromosome != hsp.sseqid:
                continue

            region_strand = normalize_strand(
                getattr(region, "strand", None)
            )

            if (
                region_strand is not None
                and hsp_strand is not None
                and region_strand != hsp_strand
            ):
                continue

            if (
                hsp_start < region.end
                and hsp_end > region.start
            ):
                region_map[idx].append(source_line)

    return region_map


def rank_proteins(
    source_lines,
    indexed_hsps,
    ssearch_hits,
):
    """
    Group SSEARCH-supported hits by parent protein.

    Score density = sum(SSEARCH bitscores) /
                    sum(SSEARCH aligned lengths).

    Only one SSEARCH alignment per original TBLASTN
    hit is counted.
    """
    groups = defaultdict(list)

    for source_line in source_lines:
        if source_line not in ssearch_hits:
            continue

        blast_hit = indexed_hsps[source_line]
        protein_id = blast_hit.qseqid

        groups[protein_id].append(
            (blast_hit, ssearch_hits[source_line])
        )

    ranked = []

    for protein_id, pairs in groups.items():
        total_bitscore = sum(
            s.bitscore for _, s in pairs
        )
        total_length = sum(
            s.alignment_length for _, s in pairs
        )

        density = total_bitscore / total_length

        ranked.append({
            "protein_id": protein_id,
            "score_density": density,
            "total_bitscore": total_bitscore,
            "aligned_length": total_length,
            "hsp_count": len(pairs),
            "blast_hsps": [b for b, _ in pairs],
            "source_lines": [
                s.source_line for _, s in pairs
            ],
        })

    ranked.sort(
        key=lambda x: (
            -x["score_density"],
            -x["total_bitscore"],
            x["protein_id"],
        )
    )

    return ranked

def summarize_ssearch_evalues(
    ranked_protein,
    ssearch_hits,
):
    """
    Retrieve SSEARCH E-values associated with
    the original TBLASTN source lines used to
    rank one protein.

    best_evalue:
        Smallest E-value among selected hits.

    worst_evalue:
        Largest E-value among selected hits.

    These values are not a combined E-value.
    """

    matched = []

    for source_line in ranked_protein["source_lines"]:
        source_line = int(source_line)

        hit = ssearch_hits.get(source_line)

        if hit is not None:
            matched.append(hit)

    if not matched:
        return {
            "best_evalue": None,
            "worst_evalue": None,
            "hsp_evalues": [],
        }

    return {
        "best_evalue": min(
            hit.evalue for hit in matched
        ),
        "worst_evalue": max(
            hit.evalue for hit in matched
        ),
        "hsp_evalues": [
            {
                "source_line": hit.source_line,
                "evalue": hit.evalue,
                "bitscore": hit.bitscore,
                "alignment_length": hit.alignment_length,
            }
            for hit in matched
        ],
    }


def annotate_regions_with_best_proteins(
    regions,
    indexed_hsps,
    ssearch_hits,
    output_path,
):
    region_map = assign_hsps_to_regions(
        regions, indexed_hsps
    )

    annotated = 0
    unmatched = 0
    missing_ssearch = 0

    output_path.parent.mkdir(
        parents=True, exist_ok=True
    )

    with output_path.open("w") as output:
        for idx, region in enumerate(regions):
            source_lines = region_map.get(idx, [])

            if not source_lines:
                unmatched += 1
                continue

            ranked = rank_proteins(
                source_lines,
                indexed_hsps,
                ssearch_hits,
            )

            if not ranked:
                missing_ssearch += 1
                continue

            best = ranked[0]

            # Preserve existing RegionAnnotation schema.
            # These ProteinHitGroup objects still contain
            # TBLASTN-derived statistics.
            protein_groups = [
                ProteinHitGroup.from_hsps(
                    item["protein_id"],
                    item["blast_hsps"],
                )
                for item in ranked
            ]

            # Preserve behavior expected by the existing
            # downstream pipeline.
            region.query_ids = [best["protein_id"]]

            best_evalues = summarize_ssearch_evalues(
                best,
                ssearch_hits,
            )

            # Preserve statistics for every ranked protein.
            all_protein_stats = []

            for item in ranked:
                evalue_stats = summarize_ssearch_evalues(
                    item,
                    ssearch_hits,
                )

                all_protein_stats.append({
                    "protein_id": item["protein_id"],
                    "score_density": item["score_density"],
                    "total_bitscore": item["total_bitscore"],
                    "aligned_length": item["aligned_length"],
                    "hsp_count": item["hsp_count"],
                    "source_lines": item["source_lines"],
                    **evalue_stats,
                })

            # Construct complete SSEARCH metadata.
            ssearch_selection = {
                "best_protein": best["protein_id"],
                "selection_method": "ssearch_score_density",
                "score_density": best["score_density"],
                "total_bitscore": best["total_bitscore"],
                "aligned_length": best["aligned_length"],
                "hsp_count": best["hsp_count"],

            # New E-value fields.
                **best_evalues,

                "all_proteins": all_protein_stats,
            }

            # Store SSEARCH metadata directly in the model.
            annotation = RegionAnnotation(
                region=region,
                best_protein=protein_groups[0],
                all_proteins=protein_groups,
                ssearch_selection=ssearch_selection,
            )

            record = annotation.to_json_dict()
            record["selection_method"] = "ssearch_score_density"

            record = annotation.to_json_dict()

            # These are the AUTHORITATIVE SSEARCH
            # statistics used for selection.
            record["selection_method"] = (
                "ssearch_score_density"
            )

            record["ssearch_selection"] = {
                "best_protein": best["protein_id"],
                "score_density": best["score_density"],
                "total_bitscore": best["total_bitscore"],
                "aligned_length": best["aligned_length"],
                "hsp_count": best["hsp_count"],
                "all_proteins": [
                    {
                        "protein_id": item["protein_id"],
                        "score_density": item["score_density"],
                        "total_bitscore": item["total_bitscore"],
                        "aligned_length": item["aligned_length"],
                        "hsp_count": item["hsp_count"],
                        "source_lines": item["source_lines"],
                    }
                    for item in ranked
                ],
            }

            # Retrieve SSEARCH hits contributing to the
            # selected parent-protein ranking.
            selected_hits = [
                ssearch_hits[int(source_line)]
                for source_line in best["source_lines"]
                if int(source_line) in ssearch_hits
            ]

            if not selected_hits:
                raise ValueError(
                    f"No SSEARCH hits for selected protein "
                    f"{best['protein_id']} in region "
                    f"{region.region_id}"
                )

            selection = record["ssearch_selection"]

            selection["best_evalue"] = min(
                hit.evalue for hit in selected_hits
            )

            selection["worst_evalue"] = max(
                hit.evalue for hit in selected_hits
            )

            selection["hsp_evalues"] = [
                {
                    "source_line": hit.source_line,
                    "evalue": hit.evalue,
                    "bitscore": hit.bitscore,
                    "alignment_length": hit.alignment_length,
                }
                for hit in selected_hits
            ]

            output.write(
                json.dumps(record) + "\n"
            )

            annotated += 1

    print(
        f"Regions: {len(regions)} | "
        f"Annotated: {annotated} | "
        f"No TBLASTN match: {unmatched} | "
        f"No SSEARCH result: {missing_ssearch}",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--regions", type=Path, required=True
    )
    parser.add_argument(
        "--tblastn", type=Path, required=True
    )
    parser.add_argument(
        "--ssearch", type=Path, required=True
    )
    parser.add_argument(
        "--output", type=Path, required=True
    )

    args = parser.parse_args()

    regions = parse_merged_bed(args.regions)

    indexed_hsps = load_indexed_tblastn(
        args.tblastn
    )

    ssearch_hits = parse_ssearch_results(
        args.ssearch
    )

    unknown_sources = (
        set(ssearch_hits) - set(indexed_hsps)
    )

    if unknown_sources:
        raise ValueError(
            "SSEARCH references TBLASTN line numbers "
            f"not present in the filtered file: "
            f"{sorted(unknown_sources)[:10]}"
        )

    annotate_regions_with_best_proteins(
        regions,
        indexed_hsps,
        ssearch_hits,
        args.output,
    )


if __name__ == "__main__":
    main()
