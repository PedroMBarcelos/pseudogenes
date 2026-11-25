#!/usr/bin/env python3
"""Audit GFF3 conventions and integrity across NCBI bacterial assemblies.

Features spanning the origin of circular sequences are supported according to
NCBI GFF3 conventions: an interval's end can exceed the declared sequence length.

The audit tests structural consistency, not biological completeness or accuracy.
It never masks sequences or modifies input files.

Python >=3.9; standard library only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

FEATURES = (
    "gene", "CDS", "pseudogene", "exon", "tRNA", "rRNA", "tmRNA",
    "repeat_region", "sequence_feature",
)


@dataclass
class Feature:
    seqid: str
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive; may exceed sequence length on circular replicon
    strand: str
    kind: str
    attrs: dict
    line: int


def open_text(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")


def choose_file(folder, patterns):
    """Choose unique plausible file, preferring a canonical filename.

    Does not confuse CDS/RNA FASTAs with the assembly genomic FASTA.
    """
    matches = []
    for pattern in patterns:
        for path in sorted(folder.glob(pattern)):
            if path.is_file() and path not in matches:
                matches.append(path)
    return matches[0] if matches else None


def attributes(raw):
    attrs = {}
    for item in raw.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs


def flagged_pseudo(attrs):
    return (
        attrs.get("pseudo", "").lower() in {"true", "1", "yes"}
        or "pseudogene" in attrs
        or attrs.get("gene_biotype", "").lower() == "pseudogene"
        or attrs.get("biotype", "").lower() == "pseudogene"
    )


def fasta_lengths(path):
    sizes = {}
    header = None
    size = 0
    with open_text(path) as stream:
        for line in stream:
            if line.startswith(">"):
                if header is not None:
                    sizes[header] = size
                header = line[1:].split()[0]
                size = 0
            elif header is not None:
                size += len(line.strip())
    if header is not None:
        sizes[header] = size
    return sizes


def read_gff(path):
    features = []
    declared_lengths = {}
    circular = set()
    meta = {}
    malformed = []
    with open_text(path) as stream:
        for lineno, line in enumerate(stream, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region "):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        lower, upper = int(parts[2]), int(parts[3])
                        if lower == 1 and upper > 0:
                            declared_lengths[parts[1]] = upper
                    except ValueError:
                        pass
            for prefix, key in (
                ("#!processor ", "processor"),
                ("#!genome-build-accession ", "assembly_accession"),
                ("#!genome-build ", "assembly_build"),
                ("#!annotation-date ", "annotation_date"),
            ):
                if line.startswith(prefix):
                    meta[key] = line[len(prefix):].strip()
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                malformed.append((lineno, "Expected 9 tab-separated fields"))
                continue
            try:
                first, last = int(fields[3]), int(fields[4])
            except ValueError:
                malformed.append((lineno, "Non-numeric feature coordinates"))
                continue
            if first < 1 or last < first:
                malformed.append((lineno, f"Invalid feature interval: {first}-{last}"))
                continue
            attrs = attributes(fields[8])
            feature = Feature(fields[0], first - 1, last, fields[6], fields[2], attrs, lineno)
            features.append(feature)
            if feature.kind == "region":
                if attrs.get("Is_circular", "").lower() in {"true", "1", "yes"}:
                    circular.add(feature.seqid)
                if feature.start == 0 and feature.seqid not in declared_lengths:
                    declared_lengths[feature.seqid] = feature.end
    return features, declared_lengths, circular, meta, malformed


def check_coordinate(feature, lengths, circular):
    """Return (reason, valid). NCBI circular-origin features may have end>L."""
    L = lengths.get(feature.seqid)
    if L is None:
        return "unknown_length", True
    if feature.start >= L:
        return "start_beyond_sequence", False
    if feature.end <= L:
        return "normal", True
    if feature.seqid in circular and feature.end - feature.start <= L:
        return "circular_origin", True
    return "end_beyond_sequence", False


def segments(feature, lengths):
    """Convert one logical feature to 1 or 2 BED-style segments."""
    L = lengths.get(feature.seqid)
    if L is not None and feature.end > L and feature.start < L:
        return [(feature.start, L), (0, feature.end - L)]
    return [(feature.start, feature.end)]


def union(intervals):
    blocks = []
    for a, b in sorted(intervals):
        if a >= b:
            continue
        if not blocks or a > blocks[-1][1]:
            blocks.append([a, b])
        elif b > blocks[-1][1]:
            blocks[-1][1] = b
    return blocks


class CoverageIndex:
    def __init__(self, targets, strand_specific, lengths):
        groups = defaultdict(list)
        for feature in targets:
            key = (feature.seqid, feature.strand) if strand_specific else feature.seqid
            groups[key].extend(segments(feature, lengths))
        self.blocks = {k: union(v) for k, v in groups.items()}
        self.ends = {k: [b for _, b in vs] for k, vs in self.blocks.items()}
        self.strand_specific = strand_specific
        self.lengths = lengths

    def fraction(self, feature):
        key = (feature.seqid, feature.strand) if self.strand_specific else feature.seqid
        blocks = self.blocks.get(key, [])
        ends = self.ends.get(key, [])
        covered = 0
        for a, b in segments(feature, self.lengths):
            i = bisect_right(ends, a)
            while i < len(blocks) and blocks[i][0] < b:
                c, d = blocks[i]
                covered += max(0, min(b, d) - max(a, c))
                i += 1
        return covered / (feature.end - feature.start)


def coverage_counts(queries, targets, same_strand, lengths):
    index = CoverageIndex(targets, same_strand, lengths)
    full = partial = none = 0
    total = 0.0
    for item in queries:
        frac = index.fraction(item)
        total += frac
        if frac >= 1 - 1e-12:
            full += 1
        elif frac > 0:
            partial += 1
        else:
            none += 1
    return full, partial, none, round(100 * total / len(queries), 2) if queries else ""


def add_issue(issues, genome, code, feature=None, note=""):
    issues.append({
        "genome": genome,
        "issue_type": code,
        "feature_type": feature.kind if feature else "",
        "seqid": feature.seqid if feature else "",
        "start_1based": feature.start + 1 if feature else "",
        "end_1based": feature.end if feature else "",
        "strand": feature.strand if feature else "",
        "line": feature.line if feature else "",
        "feature_id": feature.attrs.get("ID", "") if feature else "",
        "note": note,
    })


def audit(folder, check_fasta, issues):
    genome = folder.name
    row = {"genome": genome}
    gff_path = choose_file(folder, ["genomic.gff", "genomic.gff.gz", "genomic.gff3", "genomic.gff3.gz", "*_genomic.gff", "*_genomic.gff.gz", "*_genomic.gff3", "*_genomic.gff3.gz"])
    fasta_path = choose_file(folder, ["genomic.fna", "genomic.fna.gz", "genomic.fa", "genomic.fa.gz", "genomic.fasta", "genomic.fasta.gz", "*_genomic.fna", "*_genomic.fna.gz", "*_genomic.fa", "*_genomic.fa.gz"])
    row.update({"gff": str(gff_path) if gff_path else "", "fasta": str(fasta_path) if fasta_path else "",
                "gff_status": "PRESENT" if gff_path else "MISSING",
                "fasta_status": ("PRESENT" if fasta_path else "MISSING") if check_fasta else "NOT_CHECKED"})
    if gff_path is None:
        row["status"] = "NO_GFF"
        add_issue(issues, genome, "MISSING_GFF", note="No recognized genomic annotation GFF3 file")
        return row

    features, declared, circular, meta, malformed = read_gff(gff_path)
    row.update(meta)
    row["gff_seqids"] = len({f.seqid for f in features})
    row["circular_seqids"] = len(circular)
    row["invalid_feature_rows"] = len(malformed)
    for lineno, reason in malformed:
        add_issue(issues, genome, "MALFORMED_GFF_ROW", note=f"Line {lineno}: {reason}")

    sizes = {}
    if check_fasta:
        if fasta_path is not None:
            sizes = fasta_lengths(fasta_path)
        else:
            add_issue(issues, genome, "MISSING_FASTA", note="GFF-only analysis; cannot verify actual sequence lengths")

    seqids = {f.seqid for f in features}
    row["gff_seqids_missing_from_fasta"] = len(seqids - sizes.keys()) if sizes else ""
    row["fasta_seqids_missing_from_gff"] = len(sizes.keys() - seqids) if sizes else ""
    row["gff_fasta_length_mismatches"] = sum(declared[k] != v for k, v in sizes.items() if k in declared)
    lengths = dict(sizes)
    lengths.update(declared)  # GFF's declared coordinate system for origin-spanning features

    for seqid, size in sizes.items():
        if seqid in declared and declared[seqid] != size:
            add_issue(issues, genome, "SEQUENCE_LENGTH_DISAGREEMENT", note=f"{seqid}: GFF {declared[seqid]} nt; FASTA {size} nt")
    for seqid in sorted(seqids - sizes.keys()) if sizes else []:
        add_issue(issues, genome, "SEQID_NOT_IN_FASTA", note=seqid)
    for seqid in sorted(sizes.keys() - seqids) if sizes else []:
        add_issue(issues, genome, "FASTA_ID_NOT_IN_GFF", note=seqid)

    counts = Counter(f.kind for f in features)
    for feature in FEATURES:
        row[feature] = counts[feature]
    row["pseudo_flagged_CDS"] = sum(f.kind == "CDS" and flagged_pseudo(f.attrs) for f in features)
    row["pseudo_flagged_gene"] = sum(f.kind == "gene" and flagged_pseudo(f.attrs) for f in features)

    invalid_bounds = circular_origin = unknown_lengths = 0
    valid_features = defaultdict(list)
    for item in features:
        category, valid = check_coordinate(item, lengths, circular)
        if category == "circular_origin":
            circular_origin += 1
            add_issue(issues, genome, "CIRCULAR_ORIGIN_VALID", item,
                      note=f"end>{lengths[item.seqid]}, acceptable on a circular replicon")
        elif category == "unknown_length":
            unknown_lengths += 1
        elif not valid:
            invalid_bounds += 1
            add_issue(issues, genome, "INVALID_BOUNDS", item,
                      note=f"Sequence length={lengths[item.seqid]}, circular={item.seqid in circular}")
        if valid:
            valid_features[item.kind].append(item)

    row["circular_origin_features"] = circular_origin
    row["invalid_bounds_features"] = invalid_bounds
    row["unknown_length_features"] = unknown_lengths
    row["features_beyond_fasta_length_invalid"] = sum(
        1 for f in features
        if f.seqid in sizes and f.end > sizes[f.seqid]
        and check_coordinate(f, sizes, circular)[1] is False
    )

    gene = valid_features["gene"]
    cds = valid_features["CDS"]
    pseudo = valid_features["pseudogene"]
    rna = valid_features["tRNA"] + valid_features["rRNA"] + valid_features["tmRNA"]
    protein_genes = [f for f in gene if f.attrs.get("gene_biotype", f.attrs.get("biotype")) == "protein_coding" and not flagged_pseudo(f.attrs)]
    row["protein_coding_genes"] = len(protein_genes)
    checks = {
        "gene_by_CDS_same": (gene, cds, True),
        "protein_gene_by_CDS_same": (protein_genes, cds, True),
        "gene_by_CDS_or_RNA_same": (gene, cds + rna, True),
        "CDS_by_gene_same": (cds, gene, True),
        "CDS_by_gene_any": (cds, gene, False),
        "CDS_by_protein_gene_any": (cds, protein_genes, False),
        "RNA_by_gene_same": (rna, gene, True),
        "exon_by_RNA_same": (valid_features["exon"], rna, True),
        "pseudo_by_CDS_same": (pseudo, cds, True),
        "pseudo_by_CDS_any": (pseudo, cds, False),
        "pseudo_by_gene_any": (pseudo, gene, False),
        "pseudo_by_protein_gene_any": (pseudo, protein_genes, False),
    }
    for label, (queries, targets, same_strand) in checks.items():
        full, partial, none, mean = coverage_counts(queries, targets, same_strand, lengths)
        for suffix, value in zip(("full", "partial", "none", "mean_pct"), (full, partial, none, mean)):
            row[f"{label}_{suffix}"] = value

    fasta_mismatch = bool(sizes) and bool(
        row["gff_seqids_missing_from_fasta"] or row["fasta_seqids_missing_from_gff"]
        or row["gff_fasta_length_mismatches"] or row["features_beyond_fasta_length_invalid"]
    )
    if malformed or invalid_bounds:
        row["status"] = "FAIL_GFF_COORDINATES"
    elif fasta_mismatch:
        row["status"] = "FAIL_FASTA_MISMATCH"
    elif check_fasta and fasta_path is None:
        row["status"] = "GFF_ONLY"
    elif check_fasta and not sizes:
        row["status"] = "FAIL_EMPTY_FASTA"
    elif not check_fasta:
        row["status"] = "GFF_ONLY_UNVERIFIED"
    else:
        row["status"] = "PASS"
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("annotation_audit_156_v2.tsv"))
    parser.add_argument("--issues", type=Path, help="Optional issue and valid origin-spanning feature TSV; defaults next to --output")
    parser.add_argument("--check-fasta", action="store_true")
    args = parser.parse_args()
    if not args.root.is_dir():
        parser.error(f"No directory: {args.root}")
    folders = sorted(x for x in args.root.iterdir() if x.is_dir() and x.name.startswith(("GCA_", "GCF_")))
    if not folders:
        parser.error(f"No GCA_/GCF_ directories found in {args.root}")

    issues = []
    rows = []
    for i, folder in enumerate(folders, 1):
        try:
            row = audit(folder, args.check_fasta, issues)
        except Exception as error:
            row = {"genome": folder.name, "status": f"ERROR_{type(error).__name__}", "error": str(error)}
            add_issue(issues, folder.name, "PROCESSING_ERROR", note=f"{type(error).__name__}: {error}")
        rows.append(row)
        print(f"[{i}/{len(folders)}] {folder.name}: {row['status']}", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    headers = ["genome", "status", "gff_status", "fasta_status", "gff", "fasta"]
    headers += sorted({k for row in rows for k in row if k not in headers})
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    issue_path = args.issues or args.output.with_name(args.output.stem + "_issues.tsv")
    issue_path.parent.mkdir(parents=True, exist_ok=True)
    with issue_path.open("w", newline="") as stream:
        headers2 = ["genome", "issue_type", "feature_type", "seqid", "start_1based", "end_1based", "strand", "line", "feature_id", "note"]
        writer = csv.DictWriter(stream, fieldnames=headers2, delimiter="\t")
        writer.writeheader()
        writer.writerows(issues)

    summary = Counter(row["status"] for row in rows)
    print(f"\nAudited: {len(rows)} assembly directories")
    for status, count in sorted(summary.items()):
        print(f"  {status}: {count}")
    with_gff = [r for r in rows if r.get("gff_status") == "PRESENT"]
    with_pseudo = [r for r in with_gff if int(r.get("pseudogene") or 0) > 0]
    print(f"  GFF present: {len(with_gff)}/{len(rows)}")
    print(f"  GFFs with explicit pseudogene features: {len(with_pseudo)}/{len(with_gff)}")
    print(f"  Circular origin-spanning features (valid): {sum(int(r.get('circular_origin_features') or 0) for r in rows)}")
    print(f"  Invalid coordinate features: {sum(int(r.get('invalid_bounds_features') or 0) for r in rows)}")
    print(f"Full report: {args.output}")
    print(f"Detailed issues and circular features: {issue_path}")
    print("Reminder: PASS means structurally consistent, NOT biologically complete.")


if __name__ == "__main__":
    main()
