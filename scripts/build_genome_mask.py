#!/usr/bin/env python3
"""Build a circular-coordinate-aware 0-based BED mask from an NCBI GFF3.

Modes:
  none:                                empty BED (unmasked control)
  cds:                                 all CDS features (legacy mask)
  protein_gene:                        protein-coding gene features only
  protein_gene_preserve_pseudogene:    protein-coding genes minus ALL
                                       annotated pseudogene intervals

Pseudogene intervals include explicit 'pseudogene' features and records with
'pseudo=true', a 'pseudogene' qualifier, or a pseudogene biotype.

The program produces BED intervals; apply them with bedtools maskfasta. Do not
use the bedtools '-soft' option if hard masking is intended. The output is a
masking plan, not a measure of annotation accuracy or biological functionality.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

MODES = (
    "none", "cds", "protein_gene", "protein_gene_preserve_pseudogene"
)


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths = {}
    seqid = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                seqid = line[1:].split()[0]
                if not seqid:
                    raise ValueError("Empty sequence identifier in FASTA")
                if seqid in lengths:
                    raise ValueError(f"Duplicate FASTA identifier: {seqid}")
                lengths[seqid] = 0
            elif line.strip() and seqid is not None:
                lengths[seqid] += len(line.strip())
            elif line.strip():
                raise ValueError("FASTA sequence before first header")
    if not lengths or any(n <= 0 for n in lengths.values()):
        raise ValueError(f"Missing or empty sequences in {path}")
    return lengths


def attributes(raw: str) -> dict[str, str]:
    out = {}
    for item in raw.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
    return out


def is_pseudogene(feature_type: str, attrs: dict[str, str]) -> bool:
    biotype = attrs.get("gene_biotype", attrs.get("biotype", "")).lower()
    return (
        feature_type == "pseudogene"
        or attrs.get("pseudo", "").lower() in ("true", "1", "yes")
        or "pseudogene" in attrs
        or biotype == "pseudogene"
    )


def is_protein_gene(feature_type: str, attrs: dict[str, str]) -> bool:
    return (
        feature_type == "gene"
        and attrs.get("gene_biotype", attrs.get("biotype", "")).lower()
            == "protein_coding"
        and not is_pseudogene(feature_type, attrs)
    )


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result = []
    for start, end in sorted(intervals):
        if start >= end:
            continue
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return result


def subtract(
    intervals: list[tuple[int, int]],
    excluded: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Subtract unioned half-open excluded intervals from unioned intervals."""
    remaining = []
    j = 0
    for start, end in intervals:
        current = start
        while j < len(excluded) and excluded[j][1] <= start:
            j += 1
        k = j
        while k < len(excluded) and excluded[k][0] < end:
            ex_start, ex_end = excluded[k]
            if ex_start > current:
                remaining.append((current, min(ex_start, end)))
            current = max(current, ex_end)
            if current >= end:
                break
            k += 1
        if current < end:
            remaining.append((current, end))
    return remaining


def split_gff_feature(
    seqid: str, start: int, end: int, length: int,
    circular: bool, line_number: int
) -> list[tuple[int, int]]:
    """Convert one 1-based inclusive feature into 0-based BED parts."""
    if start < 1 or end < start or start > length:
        raise ValueError(
            f"GFF line {line_number}: invalid range {seqid}:{start}-{end} "
            f"for sequence length {length}"
        )
    if end <= length:
        return [(start - 1, end)]
    if not circular or end - start + 1 > length:
        raise ValueError(
            f"GFF line {line_number}: range {seqid}:{start}-{end} "
            f"extends beyond length {length}, but is not a valid "
            "single-turn circular-origin feature"
        )
    return [(start - 1, length), (0, end - length)]


def read_gff(path: Path, lengths: dict[str, int]):
    """Collect features and circularity before normalizing their coordinates."""
    circular = set()
    sequence_regions = {}
    records = []
    counts = Counter()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region "):
                parts = line.split()
                if len(parts) >= 4:
                    sequence_regions[parts[1]] = (int(parts[2]), int(parts[3]))
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF line {line_number}: expected 9 columns")
            seqid, _source, ftype, start, end, _score, strand, _phase, attr = fields
            try:
                start1, end1 = int(start), int(end)
            except ValueError as exc:
                raise ValueError(f"GFF line {line_number}: nonnumeric coordinates") from exc
            attrs = attributes(attr)
            if ftype == "region" and attrs.get("Is_circular", "").lower() in ("true", "1"):
                circular.add(seqid)
            records.append((line_number, seqid, ftype, start1, end1, attrs))
            counts[ftype] += 1
    groups = {"cds": defaultdict(list),
              "protein_gene": defaultdict(list),
              "pseudogene": defaultdict(list)}
    wrap_count = Counter()
    for line_number, seqid, ftype, start, end, attrs in records:
        # The GFF may list noncoding features; all still must map to the FASTA.
        if seqid not in lengths:
            raise ValueError(f"GFF line {line_number}: {seqid} is absent from FASTA")
        n = lengths[seqid]
        if seqid in sequence_regions and sequence_regions[seqid] != (1, n):
            raise ValueError(
                f"GFF ##sequence-region for {seqid} is "
                f"{sequence_regions[seqid]}, FASTA length is {n}"
            )
        parts = split_gff_feature(seqid, start, end, n, seqid in circular, line_number)
        if end > n:
            wrap_count[ftype] += 1
        if ftype == "CDS":
            groups["cds"][seqid].extend(parts)
        if is_protein_gene(ftype, attrs):
            groups["protein_gene"][seqid].extend(parts)
        if is_pseudogene(ftype, attrs):
            groups["pseudogene"][seqid].extend(parts)
    return groups, counts, wrap_count, circular


def total_bases(groups: dict[str, list[tuple[int, int]]]) -> int:
    return sum(end - start for intervals in groups.values() for start, end in intervals)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--gff", type=Path,
                        help="Required for all modes except 'none'")
    parser.add_argument("--mode", required=True, choices=MODES)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    if not args.fasta.is_file():
        parser.error(f"FASTA not found: {args.fasta}")
    if args.mode != "none" and (args.gff is None or not args.gff.is_file()):
        parser.error("This masking mode requires an existing --gff file")

    lengths = fasta_lengths(args.fasta)
    counts = Counter()
    wrapped = Counter()
    circular = set()
    selected = defaultdict(list)
    protected = defaultdict(list)
    if args.mode != "none":
        groups, counts, wrapped, circular = read_gff(args.gff, lengths)
        selected_name = "cds" if args.mode == "cds" else "protein_gene"
        selected = {chrom: merge(intervals) for chrom, intervals in
                    groups[selected_name].items()}
        protected = {chrom: merge(intervals) for chrom, intervals in
                     groups["pseudogene"].items()}
        if selected_name == "protein_gene" and not any(selected.values()):
            raise ValueError(
                "No gene features with gene_biotype=protein_coding were found. "
                "Inspect annotation conventions rather than silently producing an empty mask."
            )
    before = total_bases(selected)
    if args.mode == "protein_gene_preserve_pseudogene":
        selected = {
            chrom: subtract(intervals, protected.get(chrom, []))
            for chrom, intervals in selected.items()
        }
    selected = {chrom: merge(intervals) for chrom, intervals in selected.items()}
    after = total_bases(selected)
    args.bed.parent.mkdir(parents=True, exist_ok=True)
    with args.bed.open("w") as handle:
        for chrom in lengths:
            for start, end in selected.get(chrom, []):
                handle.write(f"{chrom}\t{start}\t{end}\n")
    summary = {
        "mode": args.mode,
        "fasta": str(args.fasta),
        "gff": str(args.gff) if args.gff else None,
        "sequences": len(lengths),
        "genome_bases": sum(lengths.values()),
        "feature_counts": dict(sorted(counts.items())),
        "circular_sequences": sorted(circular),
        "origin_spanning_features": dict(sorted(wrapped.items())),
        "mask_bases_before_pseudogene_protection": before,
        "pseudogene_protected_bases_removed_from_mask": before - after,
        "final_mask_bases": after,
        "final_mask_fraction": round(after / sum(lengths.values()), 8),
        "bed_intervals": sum(len(v) for v in selected.values()),
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
