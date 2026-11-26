#!/usr/bin/env python3

"""
Export protein-supported genomic regions to one TSV,
with three independent classification dimensions.

Input:
    JSONL from PSEUDOGENE_DETECTION

Output:
    candidate_regions.tsv

No coverage, length, identity, or mutation-based
filter is applied by this script.

Coordinates:
    BED convention, 0-based and end-exclusive.

The ORF and disruption categories are preliminary
descriptive annotations.
"""

import argparse
import csv
import json
import logging

from collections import Counter
from pathlib import Path

logger = logging.getLogger(__name__)

SHORT_THRESHOLD_NT = 300
MIN_PARTIAL_SUPPORT_NT = 30

STOP_CODONS = {"TAA", "TAG", "TGA"}

COLUMNS = [
    "region_id",
    "chromosome",
    "start",
    "end",
    "strand",
    "sequence_nt",
    "region_length_nt",
    "length_class",
    "region_type",
    "disruption_status",
    "protein_id",
    "reference_length_aa",
    "reference_coverage_pct",
    "num_hsps",
    "best_tblastn_identity_pct",
    "tblastn_total_raw_score",
    "ssearch_total_bitscore",
    "ssearch_score_density",
    "ssearch_evalue",
    "num_complete_orfs",
    "num_supported_complete_orfs",
    "best_orf_length_nt",
    "best_orf_start",
    "best_orf_end",
    "best_orf_length_class",
    "has_boundary_partial_evidence",
    "amino_acid_mismatches",
    "protein_gap_events",
    "frameshifts",
    "premature_stops",
    "missing_start_codon",
    "missing_stop_codon",
    "observed_tblastn_stop_symbols",
    "mutation_evidence",
    "tfasty_assessment_status",
    "tfasty_alignment_parsed",
    "tfasty_reference_coverage_pct",
    "tfasty_full_reference",
    "tfasty_identity_pct",
    "tfasty_similarity_pct",
    "tfasty_bitscore",
    "tfasty_evalue",
    "tfasty_query_start",
    "tfasty_query_end",
    "tfasty_alignment_code",
    "m10_alignment_count",
    "m8CC_alignment_count",
    "tfasty_failure_stage",
    "tfasty_error",
    "boundary_status",
    "boundary_error",
    "alternative_window_nt",
    "expected_start_status",
    "expected_stop_status",
    "expected_start_altered",
    "expected_stop_altered",
    "expected_start_codon",
    "expected_stop_codon",
    "expected_start_outside",
    "expected_stop_outside",
    "alternative_start_count",
    "alternative_stop_count",
]


def load_jsonl(path):
    records = []

    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue

            try:
                record = json.loads(line)

                if not isinstance(record, dict):
                    raise ValueError(
                        "Expected a JSON object"
                    )

                records.append(record)

            except (ValueError, TypeError) as exc:
                raise ValueError(
                    f"{path}:{line_number}: {exc}"
                ) from exc

    return records


def union_length_1based(intervals):
    """Union of inclusive coordinate intervals."""

    intervals = sorted(
        (min(a, b), max(a, b))
        for a, b in intervals
    )

    if not intervals:
        return 0

    start, end = intervals[0]
    total = 0

    for next_start, next_end in intervals[1:]:
        if next_start <= end + 1:
            end = max(end, next_end)
        else:
            total += end - start + 1
            start, end = next_start, next_end

    return total + end - start + 1


def reference_coverage(hsps):
    """
    Calculate reference coverage from original
    TBLASTN query intervals.

    Returns a fraction from 0 to 1.
    """

    if not hsps:
        return None

    ref_len = int(hsps[0]["qlen"])

    if ref_len <= 0:
        return None

    intervals = []

    for hsp in hsps:
        if int(hsp["qlen"]) != ref_len:
            raise ValueError(
                "Inconsistent protein lengths"
            )

        a = int(hsp["qstart"])
        b = int(hsp["qend"])

        if min(a, b) < 1 or max(a, b) > ref_len:
            raise ValueError(
                "Invalid reference alignment coordinates"
            )

        intervals.append((a, b))

    return union_length_1based(intervals) / ref_len


def length_class(length_nt):
    if length_nt is None:
        return "unknown"

    return (
        "short"
        if length_nt < SHORT_THRESHOLD_NT
        else "long"
    )


def get_supported_orfs(record):
    """
    Existing ORF detection marks protein_support
    for complete start-to-stop ORFs.
    """

    orfs = record.get("orfs") or []

    supported = [
        orf
        for orf in orfs
        if orf.get("protein_support") is True
    ]

    supported.sort(
        key=lambda orf: (
            orf.get("support_nt") or 0,
            orf.get("length_nt") or 0,
        ),
        reverse=True
    )

    return orfs, supported


def oriented_hsp_interval(hsp, region):
    """
    Convert subject coordinates to positions
    relative to the strand-oriented region.

    TBLASTN: 1-based inclusive.
    Region: BED 0-based end-exclusive.
    """

    if hsp.get("sseqid") != region["chromosome"]:
        return None

    strand = {
        "plus": "+",
        "minus": "-",
        "+": "+",
        "-": "-",
    }.get(str(hsp.get("sstrand")).lower())

    if strand != region["strand"]:
        return None

    left = min(
        int(hsp["sstart"]),
        int(hsp["send"])
    ) - 1

    right = max(
        int(hsp["sstart"]),
        int(hsp["send"])
    )

    if strand == "+":
        start = left - region["start"]
        end = right - region["start"]
    else:
        start = region["end"] - right
        end = region["end"] - left

    return start, end


def has_partial_coding_evidence(record):
    """
    Preliminary evidence for a boundary-truncated
    coding fragment supported by a TBLASTN HSP.

    Requires:
      - An extracted nucleotide sequence
      - A same-strand selected-parent alignment
      - At least 30 nt of aligned sequence
      - Alignment touching one candidate boundary
      - No in-frame stop in the aligned segment

    This does not prove a biological gene boundary
    or an actual truncation.
    """

    sequence = record.get("sequence_nt")

    if not sequence:
        return False

    sequence = sequence.upper()

    region = record["region"]
    hsps = record["best_protein"].get("hsps", [])

    for hsp in hsps:

        interval = oriented_hsp_interval(
            hsp, region
        )

        if interval is None:
            continue

        start, end = interval

        if start < 0 or end > len(sequence):
            continue

        if end - start < MIN_PARTIAL_SUPPORT_NT:
            continue

        touches_boundary = (
            start <= 3
            or len(sequence) - end <= 3
        )

        if not touches_boundary:
            continue

        aligned_dna = sequence[start:end]

        codons = [
            aligned_dna[i:i + 3]
            for i in range(
                0, len(aligned_dna) - 2, 3
            )
        ]

        if not codons:
            continue

        if any(codon in STOP_CODONS for codon in codons):
            continue

        return True

    return False


def classify_structure(record):
    """
    Conservative structural classification.

    Complete ORF requires a start-to-stop ORF
    supported by the parent-protein alignment.

    Without that support, do not infer a partial
    biological ORF from alignment boundaries.
    """

    sequence = record.get("sequence_nt")

    if not sequence:
        return "unresolved", False

    orfs, supported = get_supported_orfs(record)

    if supported:
        return "complete_orf", False

    return "unresolved", False

def classify_disruption(record):
    """
    Distinguish assessed disruption counts from
    preliminary unvalidated observations.

    None is NOT equivalent to zero.
    """

    fields = [
        "frameshifts",
        "premature_stop_codons",
        "missing_start_codon",
        "missing_stop_codon",
    ]

    values = [
        record.get(field)
        for field in fields
    ]

    # Only positive values from the validated
    # disablement fields count as possible evidence.
    if any(
        value is not None and value > 0
        for value in values
    ):
        return "possible_disruption"

    evidence = record.get("mutation_evidence") or {}

    observed_stops = (
        evidence.get("observed_tblastn_stop_symbols")
        or 0
    )

    upstream_stops = (
        evidence.get("upstream_internal_stops")
        or 0
    )

    preliminary_frameshifts = (
        evidence.get("possible_frameshifts")
    )

    # Preliminary evidence needs validation.
    if (
        observed_stops > 0
        or upstream_stops > 0
        or (
            preliminary_frameshifts is not None
            and preliminary_frameshifts > 0
        )
    ):
        return "ambiguous"

    # A genuine negative result requires all
    # relevant assessments to have been completed.
    if all(
        value is not None and value == 0
        for value in values
    ):
        return "no_disruption_detected"

    return "not_assessed"


def tfasty_columns(mutation):
    t = mutation.get("tfasty") or {}
    m = t.get("alignment_metadata") or {}
    b = mutation.get("terminal_boundary") or {}
    bs, bt = b.get("start") or {}, b.get("stop") or {}
    def pct(key):
        return 100 * float(m[key]) if m.get(key) is not None else None
    def site(d, key):
        return (d.get("expected_site") or {}).get(key)
    return {
        "tfasty_assessment_status": mutation.get("assessment_status") if t else None,
        "tfasty_alignment_parsed": "query_alignment" in t if t else None,
        "tfasty_reference_coverage_pct": t.get("reference_coverage_pct"),
        "tfasty_full_reference": t.get("full_reference"),
        "tfasty_identity_pct": pct("sw_ident"),
        "tfasty_similarity_pct": pct("sw_sim"),
        "tfasty_bitscore": m.get("tfy_bits"),
        "tfasty_evalue": m.get("tfy_expect"),
        "tfasty_query_start": t.get("query_start"),
        "tfasty_query_end": t.get("query_end"),
        "tfasty_alignment_code": t.get("alignment_code"),
        "m10_alignment_count": t.get("m10_alignment_count"),
        "m8CC_alignment_count": t.get("m8CC_alignment_count"),
        "tfasty_failure_stage": t.get("failure_stage"),
        "tfasty_error": t.get("error"),
        "boundary_status": b.get("status"),
        "boundary_error": b.get("reason"),
        "alternative_window_nt": b.get("alternative_window_nt"),
        "expected_start_status": bs.get("expected_site_status"),
        "expected_stop_status": bt.get("expected_site_status"),
        "expected_start_altered": bs.get("expected_site_altered"),
        "expected_stop_altered": bt.get("expected_site_altered"),
        "expected_start_codon": site(bs, "codon"),
        "expected_stop_codon": site(bt, "codon"),
        "expected_start_outside": site(bs, "outside_candidate"),
        "expected_stop_outside": site(bt, "outside_candidate"),
        "alternative_start_count": len(bs.get("alternatives", [])) if bs else None,
        "alternative_stop_count": len(bt.get("alternatives", [])) if bt else None,
    }


def record_to_row(record):
    region = record["region"]
    protein = record["best_protein"]
    hsps = protein.get("hsps", [])

    ssearch = record.get("ssearch_selection") or {}
    mutation = record.get("mutation_evidence") or {}

    region_id = region["region_id"]

    region_start = int(region["start"])
    region_end = int(region["end"])

    if region_start < 0 or region_end <= region_start:
        raise ValueError(
            f"{region_id}: invalid region coordinates"
        )

    sequence = record.get("sequence_nt")

    if sequence is not None:
        sequence = sequence.upper()

        if len(sequence) != region_end - region_start:
            raise ValueError(
                f"{region_id}: sequence length does not "
                "match genomic interval"
            )

    region_length = region_end - region_start

    structure, partial_evidence = (
        classify_structure(record)
    )

    disruption = classify_disruption(record)

    orfs, supported_orfs = get_supported_orfs(record)

    best_orf = (
        supported_orfs[0]
        if supported_orfs else None
    )

    best_orf_length = (
        best_orf.get("length_nt")
        if best_orf else None
    )

    coverage = reference_coverage(hsps)

    coverage_pct = (
        round(coverage * 100, 3)
        if coverage is not None else None
    )

    identities = [
        float(hsp["pident"])
        for hsp in hsps
    ]

    # This is the maximum individual TBLASTN HSP
    # identity, not full-region identity.
    best_identity = (
        max(identities)
        if identities else None
    )

    # Optional SSEARCH value: only present if
    # saved by the upstream annotation parser.
    ssearch_evalue = next(
        (
            ssearch[key]
            for key in (
                "best_evalue",
                "min_evalue",
                "evalue",
            )
            if ssearch.get(key) is not None
        ),
        None
    )

    return {
        "region_id": region_id,
        "chromosome": region["chromosome"],
        "start": region_start,
        "end": region_end,
        "strand": region["strand"],
        "sequence_nt": sequence,
        "region_length_nt": region_length,
        "length_class": length_class(region_length),
        "region_type": structure,
        "disruption_status": disruption,
        "protein_id": protein["protein_id"],
        "reference_length_aa":
            int(hsps[0]["qlen"]) if hsps else None,
        "reference_coverage_pct": coverage_pct,
        "num_hsps": len(hsps),
        "best_tblastn_identity_pct": best_identity,
        "tblastn_total_raw_score":
            protein.get("total_score"),
        "ssearch_total_bitscore":
            ssearch.get("total_bitscore"),
        "ssearch_score_density":
            ssearch.get("score_density"),
        "ssearch_evalue": ssearch_evalue,
        "num_complete_orfs": len(orfs),
        "num_supported_complete_orfs":
            len(supported_orfs),
        "best_orf_length_nt": best_orf_length,
        "best_orf_start":
            best_orf.get("genomic_start")
            if best_orf else None,
        "best_orf_end":
            best_orf.get("genomic_end")
            if best_orf else None,
        "best_orf_length_class":
            length_class(best_orf_length)
            if best_orf else "not_detected",
        "has_boundary_partial_evidence":
            partial_evidence,
        "amino_acid_mismatches":
            record.get("non_synonymous_substitutions"),
        "protein_gap_events":
            record.get("in_frame_indels"),
        "frameshifts":
            record.get("frameshifts"),
        "premature_stops":
            record.get("premature_stop_codons"),
        "missing_start_codon":
            record.get("missing_start_codon"),
        "missing_stop_codon":
            record.get("missing_stop_codon"),
        "observed_tblastn_stop_symbols":
            mutation.get("observed_tblastn_stop_symbols"),
        **tfasty_columns(mutation),
        "mutation_evidence": json.dumps(
            mutation,
            ensure_ascii=False
        ),
    }


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input_jsonl",
        type=Path
    )

    parser.add_argument(
        "output_tsv",
        type=Path
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s"
    )

    records = load_jsonl(args.input_jsonl)

    rows = [
        record_to_row(record)
        for record in records
    ]

    args.output_tsv.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with args.output_tsv.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=COLUMNS,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(rows)

    logger.info(
        "Exported %d regions to %s",
        len(rows),
        args.output_tsv
    )

    for column in (
        "region_type",
        "length_class",
        "disruption_status",
    ):
        logger.info(
            "%s: %s",
            column,
            dict(
                Counter(
                    row[column]
                    for row in rows
                )
            )
        )

    missing_dna = sum(
        not row["sequence_nt"]
        for row in rows
    )

    if missing_dna:
        logger.warning(
            "%d regions have no extracted DNA. "
            "Check the characterization stage.",
            missing_dna
        )


if __name__ == "__main__":
    main()
