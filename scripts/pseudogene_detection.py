
#!/usr/bin/env python3

"""
Pseudogene candidate characterization.

Pipeline:
    TBLASTN
        -> Candidate genomic regions
        -> SSEARCH-based parent protein selection
        -> Pseudogene candidate characterization

This script:
    1. Loads RegionAnnotation objects from JSONL.
    2. Preserves the SSEARCH-selected parent protein.
    3. Validates candidate coordinates against the original genome.
    4. Calculates reference protein coverage from TBLASTN HSPs.
    5. Calculates descriptive protein-alignment differences.
    6. Preserves unknown gene-disruption measurements as None.
    7. Assigns a conservative candidate classification.
    8. Writes PseudogeneAnnotation objects to JSONL.

IMPORTANT:
    - TBLASTN HSPs are used for genomic mapping and
      descriptive alignment statistics.
    - SSEARCH is used for parent protein selection.
    - Frameshifts are NOT inferred from distinct TBLASTN frames.
    - Missing start/stop codons are NOT inferred from
      incomplete alignments.
    - Protein gaps are NOT evidence of nucleotide frameshifts.
    - None means not assessed, not zero.

Expected arguments:
    python3 pseudogene_detection.py \
        region_annotations.jsonl \
        original_genome.fasta \
        initial_pseudogene_annotations.jsonl
"""

import argparse
import json
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Tuple

from Bio import SeqIO

from models.models import (
    BlastHit,
    RegionAnnotation,
    DisablementCounts,
    PseudogeneAnnotation,
)


logger = logging.getLogger(__name__)


# ============================================================
# CONFIGURATION
# ============================================================

GAP_CHAR = "-"
STOP_CHAR = "*"

# Descriptive region size threshold.
# This does not establish whether a true ORF is short.
MIN_REGION_SIZE_NT = 300

# Used only for describing partial protein alignments.
# This is NOT a pseudogene classification threshold.
#MIN_REFERENCE_COVERAGE = 0.80


# ============================================================
# LOAD REGION ANNOTATIONS
# ============================================================

def load_region_annotations_from_json(
    annotations_path: Path,
) -> List[RegionAnnotation]:
    """
    Load RegionAnnotation objects from JSONL.

    The updated RegionAnnotation.from_dict() preserves
    the SSEARCH selection information.
    """

    annotations = []

    with annotations_path.open(
        "r", encoding="utf-8"
    ) as handle:

        for line_number, line in enumerate(
            handle, start=1
        ):
            line = line.strip()

            if not line:
                continue

            try:
                record = json.loads(line)

                annotation = (
                    RegionAnnotation.from_dict(record)
                )

            except (ValueError, TypeError, KeyError) as exc:
                raise ValueError(
                    f"Cannot parse annotation at "
                    f"{annotations_path}:{line_number}"
                ) from exc

            annotations.append(annotation)

    logger.info(
        "Loaded %d region annotations from %s",
        len(annotations),
        annotations_path,
    )

    return annotations


# ============================================================
# VALIDATE GENOMIC REGION
# ============================================================

def validate_genomic_region(
    region,
    genome_lengths,
) -> None:
    """
    Validate a candidate against the original genome.

    Assumes BED coordinates:
        start: 0-based inclusive
        end:   0-based exclusive

    Does not infer gene boundaries or coding frames.
    """

    chromosome = region.chromosome

    if chromosome not in genome_lengths:
        raise ValueError(
            f"Region {region.region_id}: "
            f"chromosome '{chromosome}' was not found "
            f"in the original genome FASTA."
        )

    chromosome_length = genome_lengths[chromosome]

    if not (
        0 <= region.start < region.end
        <= chromosome_length
    ):
        raise ValueError(
            f"Region {region.region_id}: invalid "
            f"coordinates {region.start}-{region.end} "
            f"for chromosome {chromosome} "
            f"(length {chromosome_length})."
        )

    if region.strand not in ("+", "-"):
        raise ValueError(
            f"Region {region.region_id}: "
            f"invalid strand '{region.strand}'."
        )


def load_genome_lengths(genome_fasta_path: Path):
    """
    Read original genome FASTA sequence lengths.

    The current script needs coordinates and lengths,
    not full genome sequences.

    A later frameshift-aware stage will access
    the actual nucleotide sequences.
    """

    lengths = {}

    for record in SeqIO.parse(
        str(genome_fasta_path), "fasta"
    ):
        if record.id in lengths:
            raise ValueError(
                f"Duplicate chromosome ID: {record.id}"
            )

        lengths[record.id] = len(record.seq)

    if not lengths:
        raise ValueError(
            f"No sequences found in {genome_fasta_path}"
        )

    logger.info(
        "Loaded %d genomic sequence lengths",
        len(lengths),
    )

    return lengths


# ============================================================
# VALIDATE SELECTED PARENT PROTEIN
# ============================================================

def validate_selected_protein(
    region_annotation: RegionAnnotation,
) -> None:
    """
    Ensure the stored HSPs belong to the parent selected
    for this genomic region.

    The parent identity comes from SSEARCH selection.

    The HSPs remain original TBLASTN alignments.
    """

    best_protein = region_annotation.best_protein

    selected_id = best_protein.protein_id

    selection = region_annotation.ssearch_selection

    if selection is not None:
        ssearch_id = selection.get("best_protein")

        if (
            ssearch_id is not None
            and ssearch_id != selected_id
        ):
            raise ValueError(
                f"Region {region_annotation.region.region_id}: "
                f"SSEARCH parent '{ssearch_id}' differs "
                f"from stored parent '{selected_id}'."
            )
    else:
        logger.warning(
            "Region %s has no SSEARCH selection metadata. "
            "Using the stored best_protein.",
            region_annotation.region.region_id,
        )

    # The original region.parent_protein may have been
    # assigned before SSEARCH-based parent selection.
    initial_parent = (
        region_annotation.region.parent_protein
    )

    if (
        initial_parent
        and initial_parent != selected_id
    ):
        logger.debug(
            "Region %s: original parent %s differs "
            "from SSEARCH-selected parent %s",
            region_annotation.region.region_id,
            initial_parent,
            selected_id,
        )

    hsps = best_protein.hsps

    for hsp in hsps:
        if hsp.qseqid != selected_id:
            raise ValueError(
                f"Region {region_annotation.region.region_id}: "
                f"HSP from protein {hsp.qseqid} belongs "
                f"to selected parent {selected_id}."
            )

    if hsps:
        reference_lengths = {
            hsp.qlen for hsp in hsps
        }

        if len(reference_lengths) != 1:
            raise ValueError(
                f"Region {region_annotation.region.region_id}: "
                "inconsistent reference protein lengths."
            )


# ============================================================
# REFERENCE PROTEIN COVERAGE
# ============================================================

def calculate_reference_coverage(
    hsps: List[BlastHit],
) -> Optional[float]:
    """
    Calculate reference protein coverage.

    Coverage is the union of aligned QUERY intervals
    divided by reference protein length.

    TBLASTN query coordinates:
        1-based inclusive

    This avoids double-counting overlapping HSPs.

    Returns:
        Coverage as a fraction from 0 to 1.

        None if coverage cannot be calculated.
    """

    if not hsps:
        return None

    reference_length = hsps[0].qlen

    if reference_length <= 0:
        return None

    intervals = []

    for hsp in hsps:

        if hsp.qlen != reference_length:
            raise ValueError(
                "Inconsistent reference protein lengths."
            )

        start = min(hsp.qstart, hsp.qend)
        end = max(hsp.qstart, hsp.qend)

        if not (
            1 <= start <= end <= reference_length
        ):
            raise ValueError(
                f"Invalid query interval: "
                f"{start}-{end}, qlen={reference_length}"
            )

        intervals.append((start, end))

    intervals.sort()

    merged = []

    for start, end in intervals:

        if not merged:
            merged.append([start, end])
            continue

        previous = merged[-1]

        if start <= previous[1] + 1:
            previous[1] = max(
                previous[1], end
            )
        else:
            merged.append([start, end])

    covered_residues = sum(
        end - start + 1
        for start, end in merged
    )

    coverage = covered_residues / reference_length

    return coverage


# ============================================================
# AMINO-ACID MISMATCHES
# ============================================================

def count_amino_acid_mismatches(
    hsps: List[BlastHit],
) -> int:
    """
    Count amino-acid mismatches in TBLASTN HSPs.

    Excludes:
        - Gap positions
        - Stop symbols
        - Unknown residues (X)

    IMPORTANT:
        This is a descriptive alignment measurement.

        It is not a count of confirmed gene-disabling
        nonsynonymous nucleotide mutations.

        Overlapping HSPs can contribute duplicate
        alignment evidence.
    """

    total_mismatches = 0

    for hsp in hsps:

        query_seq = hsp.qseq.upper()
        subject_seq = hsp.sseq.upper()

        if len(query_seq) != len(subject_seq):
            logger.warning(
                "Unequal aligned sequence lengths for "
                "protein %s: query=%d, subject=%d. "
                "Skipping this HSP's mismatch count.",
                hsp.qseqid,
                len(query_seq),
                len(subject_seq),
            )
            continue

        for q_aa, s_aa in zip(
            query_seq, subject_seq
        ):

            if (
                GAP_CHAR in (q_aa, s_aa)
                or STOP_CHAR in (q_aa, s_aa)
                or "X" in (q_aa, s_aa)
            ):
                continue

            if q_aa != s_aa:
                total_mismatches += 1

    return total_mismatches


# ============================================================
# PROTEIN ALIGNMENT GAP EVENTS
# ============================================================

def count_protein_gap_events(
    hsps: List[BlastHit],
) -> int:
    """
    Count contiguous gap runs in protein alignments.

    Example:
        "AA---BB" contains one gap run.

    These gaps describe amino-acid-level alignments.

    They are NOT verified nucleotide indels and
    must NOT be interpreted as frameshifts.
    """

    total_events = 0

    for hsp in hsps:

        query_gap_runs = re.findall(
            r"-+", hsp.qseq
        )

        subject_gap_runs = re.findall(
            r"-+", hsp.sseq
        )

        total_events += (
            len(query_gap_runs)
            + len(subject_gap_runs)
        )

    return total_events


# ============================================================
# DESCRIPTIVE ALIGNMENT STATISTICS
# ============================================================

def calculate_alignment_statistics(
    hsps: List[BlastHit],
):
    """
    Return descriptive statistics for the selected
    parent protein's TBLASTN HSPs.

    These measurements are not necessarily independent
    because HSPs may overlap.
    """

    return {
        "amino_acid_mismatches":
            count_amino_acid_mismatches(hsps),

        "protein_gap_events":
            count_protein_gap_events(hsps),

        "reference_coverage":
            calculate_reference_coverage(hsps),

        "num_hsps":
            len(hsps),
    }


# ============================================================
# DISABLEMENT ASSESSMENT
# ============================================================

def analyze_protein_for_disablements(
    region_annotation: RegionAnnotation,
) -> Tuple[DisablementCounts, Optional[float]]:
    """
    Characterize the selected protein alignment.

    Keep unassessed genomic disruptions as None.

    The present pipeline has not yet established
    reliable nucleotide-level evidence for:

        - Frameshifts
        - Missing initiation codons
        - Missing termination codons
        - Confirmed internal stop codons
        - Gene-length disruption

    The region-level internal_stops field, if present,
    remains available in region_annotation.region,
    but is not automatically promoted to a verified
    mutation in DisablementCounts.

    This avoids treating HSP-level observations from
    competing proteins as validated gene disruptions.
    """

    hsps = region_annotation.best_protein.hsps

    stats = calculate_alignment_statistics(hsps)

    counts = DisablementCounts(

        # Legacy field names preserved for downstream
        # serialization compatibility.
        #
        # These are descriptive alignment statistics.
        non_synonymous_substitutions=(
            stats["amino_acid_mismatches"]
        ),

        in_frame_indels=(
            stats["protein_gap_events"]
        ),

        # The following measurements have not been
        # reliably assessed yet.
        frameshifts=None,

        missing_start_codon=None,

        missing_stop_codon=None,

        premature_stop_codons=None,

        size_mismatch=None,
    )

    return (
        counts,
        stats["reference_coverage"],
    )


# ============================================================
# CLASSIFICATION
# ============================================================

def classify_candidate(
    disablements: DisablementCounts,
) -> Tuple[Optional[bool], str]:
    """
    Classify candidate using validated evidence only.

    Returns:
        (True, "pseudogene")
            At least one verified disabling mutation.

        (False, "no_disruption_detected")
            All required checks completed,
            none detected.

        (None, "unresolved")
            Evidence is incomplete.

    NOTE:
        This function assumes any non-None disruption
        measurements have already been validated.

        Do not populate these fields from TBLASTN
        frame counts or protein alignment gaps.
    """

    evidence = (
        disablements.frameshifts,
        disablements.premature_stop_codons,
        disablements.missing_start_codon,
        disablements.missing_stop_codon,
        disablements.size_mismatch,
    )

    # At least one verified disabling event.
    if any(
        value is not None and value > 0
        for value in evidence
    ):
        return True, "pseudogene"

    # At least one assessment is missing.
    if any(value is None for value in evidence):
        return None, "unresolved"

    # Every required measurement was performed.
    # No disabling event was detected.
    return False, "no_disruption_detected"


# ============================================================
# PROCESS SINGLE REGION
# ============================================================


# ============================================================
# ORF DETECTION
# ============================================================

SMALL_ORF_MAX_NT = 300
MIN_ORF_SUPPORT_NT = 30

# Starting assumptions for bacterial genomes.
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TGA", "TAG"}


def find_complete_orfs(sequence_nt: str):
    """
    Find complete start-to-stop ORFs in the
    strand-oriented nucleotide sequence.

    Returns 0-based, end-exclusive coordinates
    relative to the extracted sequence.

    Includes the terminal stop codon in ORF length.
    """

    sequence_nt = sequence_nt.upper()
    orfs = []

    for frame in range(3):
        starts = []

        for pos in range(
            frame,
            len(sequence_nt) - 2,
            3
        ):
            codon = sequence_nt[pos:pos + 3]

            if codon in STOP_CODONS:
                for start in starts:
                    orfs.append({
                        "start": start,
                        "end": pos + 3,
                        "frame": frame,
                        "length_nt": pos + 3 - start,
                        "protein_support": False,
                        "support_nt": 0,
                        "supporting_hsps": 0,
                    })

                starts = []

            elif codon in START_CODONS:
                starts.append(pos)

    return orfs


def hsp_oriented_interval(hsp, region):
    """
    Convert a TBLASTN genomic HSP interval to
    coordinates relative to the strand-oriented
    candidate sequence.

    Region:
        BED coordinates (0-based, end-exclusive).

    HSP:
        TBLASTN coordinates (1-based, inclusive).
    """

    if hsp.sseqid != region.chromosome:
        return None

    strand = {
        "plus": "+",
        "minus": "-",
        "+": "+",
        "-": "-",
    }.get(str(hsp.sstrand).lower())

    if strand != region.strand:
        return None

    genomic_start = min(hsp.sstart, hsp.send) - 1
    genomic_end = max(hsp.sstart, hsp.send)

    if region.strand == "+":
        start = genomic_start - region.start
        end = genomic_end - region.start

    else:
        start = region.end - genomic_end
        end = region.end - genomic_start

    # HSP must be fully represented by the region.
    if start < 0 or end > region.length:
        return None

    return start, end


def interval_union_length(intervals):
    """Total unique nucleotide positions covered."""

    if not intervals:
        return 0

    intervals = sorted(intervals)

    total = 0
    start, end = intervals[0]

    for next_start, next_end in intervals[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end

    total += end - start
    return total


def evaluate_orf_protein_support(
    orfs,
    region_annotation,
):
    """
    Evaluate positional and reading-frame consistency
    between complete ORFs and selected-parent HSPs.

    This is preliminary protein support, not a new
    independent sequence-homology validation.
    """

    region = region_annotation.region
    hsps = region_annotation.best_protein.hsps

    hsp_intervals = []

    for hsp in hsps:
        interval = hsp_oriented_interval(hsp, region)

        if interval is not None:
            hsp_intervals.append(interval)

    for orf in orfs:

        # Exclude the terminal stop codon from
        # the protein-coding interval.
        coding_end = orf["end"] - 3
        supported_intervals = []
        supporting_hsps = 0

        for hsp_start, hsp_end in hsp_intervals:

            overlap_start = max(
                orf["start"], hsp_start
            )
            overlap_end = min(
                coding_end, hsp_end
            )

            overlap_nt = max(
                0, overlap_end - overlap_start
            )

            hsp_length = hsp_end - hsp_start

            # Require sufficient overlap with HSP.
            if overlap_nt < MIN_ORF_SUPPORT_NT:
                continue

            if overlap_nt < 0.5 * hsp_length:
                continue

            # ORF and HSP must share a reading frame.
            if (
                hsp_start - orf["start"]
            ) % 3 != 0:
                continue

            supported_intervals.append(
                (overlap_start, overlap_end)
            )

            supporting_hsps += 1

        support_nt = interval_union_length(
            supported_intervals
        )

        orf["support_nt"] = support_nt
        orf["supporting_hsps"] = supporting_hsps
        orf["protein_support"] = support_nt > 0

    return orfs



# ============================================================
# ORF DETECTION
# ============================================================

SMALL_ORF_MAX_NT = 300
MIN_ORF_SUPPORT_NT = 30

# Starting assumptions for bacterial genomes.
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TGA", "TAG"}


def find_complete_orfs(sequence_nt: str):
    """
    Find complete start-to-stop ORFs in the
    strand-oriented nucleotide sequence.

    Returns 0-based, end-exclusive coordinates
    relative to the extracted sequence.

    Includes the terminal stop codon in ORF length.
    """

    sequence_nt = sequence_nt.upper()
    orfs = []

    for frame in range(3):
        starts = []

        for pos in range(
            frame,
            len(sequence_nt) - 2,
            3
        ):
            codon = sequence_nt[pos:pos + 3]

            if codon in STOP_CODONS:
                for start in starts:
                    orfs.append({
                        "start": start,
                        "end": pos + 3,
                        "frame": frame,
                        "length_nt": pos + 3 - start,
                        "protein_support": False,
                        "support_nt": 0,
                        "supporting_hsps": 0,
                    })

                starts = []

            elif codon in START_CODONS:
                starts.append(pos)

    return orfs


def hsp_oriented_interval(hsp, region):
    """
    Convert a TBLASTN genomic HSP interval to
    coordinates relative to the strand-oriented
    candidate sequence.

    Region:
        BED coordinates (0-based, end-exclusive).

    HSP:
        TBLASTN coordinates (1-based, inclusive).
    """

    if hsp.sseqid != region.chromosome:
        return None

    strand = {
        "plus": "+",
        "minus": "-",
        "+": "+",
        "-": "-",
    }.get(str(hsp.sstrand).lower())

    if strand != region.strand:
        return None

    genomic_start = min(hsp.sstart, hsp.send) - 1
    genomic_end = max(hsp.sstart, hsp.send)

    if region.strand == "+":
        start = genomic_start - region.start
        end = genomic_end - region.start

    else:
        start = region.end - genomic_end
        end = region.end - genomic_start

    # HSP must be fully represented by the region.
    if start < 0 or end > region.length:
        return None

    return start, end


def interval_union_length(intervals):
    """Total unique nucleotide positions covered."""

    if not intervals:
        return 0

    intervals = sorted(intervals)

    total = 0
    start, end = intervals[0]

    for next_start, next_end in intervals[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end

    total += end - start
    return total


def evaluate_orf_protein_support(
    orfs,
    region_annotation,
):
    """
    Evaluate positional and reading-frame consistency
    between complete ORFs and selected-parent HSPs.

    This is preliminary protein support, not a new
    independent sequence-homology validation.
    """

    region = region_annotation.region
    hsps = region_annotation.best_protein.hsps

    hsp_intervals = []

    for hsp in hsps:
        interval = hsp_oriented_interval(hsp, region)

        if interval is not None:
            hsp_intervals.append(interval)

    for orf in orfs:

        # Exclude the terminal stop codon from
        # the protein-coding interval.
        coding_end = orf["end"] - 3
        supported_intervals = []
        supporting_hsps = 0

        for hsp_start, hsp_end in hsp_intervals:

            overlap_start = max(
                orf["start"], hsp_start
            )
            overlap_end = min(
                coding_end, hsp_end
            )

            overlap_nt = max(
                0, overlap_end - overlap_start
            )

            hsp_length = hsp_end - hsp_start

            # Require sufficient overlap with HSP.
            if overlap_nt < MIN_ORF_SUPPORT_NT:
                continue

            if overlap_nt < 0.5 * hsp_length:
                continue

            # ORF and HSP must share a reading frame.
            if (
                hsp_start - orf["start"]
            ) % 3 != 0:
                continue

            supported_intervals.append(
                (overlap_start, overlap_end)
            )

            supporting_hsps += 1

        support_nt = interval_union_length(
            supported_intervals
        )

        orf["support_nt"] = support_nt
        orf["supporting_hsps"] = supporting_hsps
        orf["protein_support"] = support_nt > 0

    return orfs



def characterize_region(
    region_annotation: RegionAnnotation,
    genome_seqs,
    genome_lengths,
) -> PseudogeneAnnotation:
    """
    Characterize a protein-associated genomic region.

    Extract DNA, measure length and reference coverage,
    detect complete ORFs, examine preliminary mutation
    evidence, and classify as small_orf or candidate.

    Coverage and mutation observations never cause
    a candidate to be discarded.
    """

    region = region_annotation.region

    # --------------------------------------------------------
    # 1. VALIDATE REGION
    # --------------------------------------------------------

    validate_genomic_region(
        region,
        genome_lengths,
    )

    validate_selected_protein(
        region_annotation
    )

    # --------------------------------------------------------
    # 2. EXTRACT ORIGINAL GENOMIC SEQUENCE
    # --------------------------------------------------------

    dna = genome_seqs[
        region.chromosome
    ][region.start:region.end]

    if region.strand == "-":
        dna = dna.reverse_complement()

    sequence_nt = str(dna).upper()
    region_length_nt = len(sequence_nt)

    if region_length_nt != region.length:
        raise ValueError(
            f"Region {region.region_id}: "
            "extracted sequence length differs "
            "from BED region length."
        )

    # --------------------------------------------------------
    # 3. CALCULATE ALIGNMENT MEASUREMENTS
    # --------------------------------------------------------

    disablements, reference_coverage = (
        analyze_protein_for_disablements(
            region_annotation
        )
    )

    # Descriptive only; not a filtering criterion.
    is_small_region = (
        region_length_nt < SMALL_ORF_MAX_NT
    )

    # --------------------------------------------------------
    # 4. FIND COMPLETE ORFS
    # --------------------------------------------------------

    orfs = find_complete_orfs(
        sequence_nt
    )

    orfs = evaluate_orf_protein_support(
        orfs,
        region_annotation
    )

    # Add original genomic coordinates to each ORF.
    for orf in orfs:

        if region.strand == "+":
            orf["genomic_start"] = (
                region.start + orf["start"]
            )
            orf["genomic_end"] = (
                region.start + orf["end"]
            )

        else:
            orf["genomic_start"] = (
                region.end - orf["end"]
            )
            orf["genomic_end"] = (
                region.end - orf["start"]
            )

    supported_orfs = [
        orf for orf in orfs
        if orf["protein_support"]
    ]

    # Choose the ORF with the strongest genomic
    # overlap with selected-parent alignments.
    # Do not prefer small ORFs automatically.
    supported_orfs.sort(
        key=lambda orf: (
            orf["support_nt"],
            orf["length_nt"]
        ),
        reverse=True
    )

    best_orf = (
        supported_orfs[0]
        if supported_orfs
        else None
    )

    best_orf_length_nt = (
        best_orf["length_nt"]
        if best_orf is not None
        else None
    )

    # --------------------------------------------------------
    # 5. TWO-CLASS CLASSIFICATION
    # --------------------------------------------------------

    if (
        best_orf is not None
        and best_orf_length_nt <= SMALL_ORF_MAX_NT
    ):
        classification_status = "small_orf"
        is_small_orf = True

    else:
        classification_status = "candidate"
        is_small_orf = False

    # We are not assigning pseudogene/functional
    # status in this version.
    is_pseudogene = None

    # --------------------------------------------------------
    # 6. RECORD POSSIBLE MUTATION EVIDENCE
    # --------------------------------------------------------

    hsps = region_annotation.best_protein.hsps

    mutation_evidence = {
        "amino_acid_mismatches":
            disablements.non_synonymous_substitutions,

        "protein_gap_events":
            disablements.in_frame_indels,

        # Observed in translated TBLASTN alignments.
        # Not necessarily independent mutation events.
        "observed_tblastn_stop_symbols": sum(
            hsp.sseq.count("*")
            for hsp in hsps
        ),

        # Retain earlier region-level evidence
        # without treating it as confirmed.
        "upstream_internal_stops":
            region.internal_stops,

        "possible_frameshifts":
            None,

        "missing_start_codon":
            None,

        "missing_stop_codon":
            None,

        "assessment_status":
            "preliminary",
    }

    # --------------------------------------------------------
    # 7. CREATE ANNOTATION
    # --------------------------------------------------------

    return PseudogeneAnnotation(
        region_annotation=region_annotation,

        disablements=disablements,

        sequence_nt=sequence_nt,
        region_length_nt=region_length_nt,

        reference_coverage=reference_coverage,

        is_small_region=is_small_region,

        orfs=orfs,
        best_orf_length_nt=best_orf_length_nt,

        mutation_evidence=mutation_evidence,

        classification_status=classification_status,
        is_small_orf=is_small_orf,
        is_pseudogene=is_pseudogene,
    )

# ============================================================
# PROCESS ALL CANDIDATE REGIONS
# ============================================================

def annotate_pseudogenes(
    region_annotations: List[RegionAnnotation],
    genome_fasta_path: Path,
    output_path: Path,
) -> List[PseudogeneAnnotation]:
    """
    Characterize all candidate genomic regions.

    Preserve unresolved candidates in the output.

    Output:
        One JSON object per candidate region.
    """

    logger.info(
        "Loading original genome information from %s",
        genome_fasta_path,
    )

    genome_seqs = {
        record.id: record.seq
        for record in SeqIO.parse(
            str(genome_fasta_path), "fasta"
        )
    }

    genome_lengths = {
        chromosome: len(sequence)
        for chromosome, sequence in genome_seqs.items()
    }

    if not genome_seqs:
        raise ValueError(
            f"No sequences found in {genome_fasta_path}"
        )    

    pseudogene_annotations = []

    for region_annotation in region_annotations:

        candidate = characterize_region(
            region_annotation,
            genome_seqs,
            genome_lengths,
        )

        pseudogene_annotations.append(candidate)

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    status_counts = defaultdict(int)

    for annotation in pseudogene_annotations:
        status_counts[
            annotation.classification_status
        ] += 1

    logger.info(
        "Total candidate regions: %d",
        len(pseudogene_annotations),
    )

    logger.info(
        "Pseudogene: %d",
        status_counts["pseudogene"],
    )

    logger.info(
        "No disruption detected: %d",
        status_counts["no_disruption_detected"],
    )

    logger.info(
        "Unresolved: %d",
        status_counts["unresolved"],
    )

    # --------------------------------------------------------
    # WRITE JSONL
    # --------------------------------------------------------

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w", encoding="utf-8"
    ) as handle:

        for annotation in pseudogene_annotations:

            record = annotation.to_json_dict()

            handle.write(
                json.dumps(
                    record,
                    ensure_ascii=False
                ) + "\n"
            )

    logger.info(
        "Saved %d candidate annotations to %s",
        len(pseudogene_annotations),
        output_path,
    )

    return pseudogene_annotations


# ============================================================
# COMMAND LINE INTERFACE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Characterize candidate pseudogene regions "
            "using SSEARCH-selected parent proteins "
            "and original TBLASTN alignment evidence."
        )
    )

    parser.add_argument(
        "region_annotations",
        type=Path,
        help="Region annotations JSONL",
    )

    parser.add_argument(
        "genome_fasta",
        type=Path,
        help="Original unmasked genome FASTA",
    )

    parser.add_argument(
        "output",
        type=Path,
        help="Output pseudogene annotations JSONL",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=(
            logging.DEBUG
            if args.verbose
            else logging.INFO
        ),
        format=(
            "%(asctime)s - %(levelname)s - %(message)s"
        ),
    )

    for input_path in (
        args.region_annotations,
        args.genome_fasta,
    ):
        if not input_path.is_file():
            parser.error(
                f"Input file not found: {input_path}"
            )

    region_annotations = (
        load_region_annotations_from_json(
            args.region_annotations
        )
    )

    annotate_pseudogenes(
        region_annotations=region_annotations,
        genome_fasta_path=args.genome_fasta,
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
