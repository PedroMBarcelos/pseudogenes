#!/usr/bin/env python3
"""Batch-audit NCBI GFF3 feature coverage for multiple bacterial assemblies.

Structural coverage is not an independent assessment of biological annotation
accuracy or pseudogene detection sensitivity.
"""
import argparse
import csv
import gzip
from bisect import bisect_right
from collections import Counter, defaultdict
from pathlib import Path

FEATURES = ("gene", "CDS", "pseudogene", "exon", "tRNA", "rRNA", "tmRNA", "repeat_region", "sequence_feature")


def opengz(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")


def get_file(folder, basename):
    for name in (basename, basename + ".gz"):
        p = folder / name
        if p.is_file():
            return p
    return None


def attr_dict(raw):
    out = {}
    for part in raw.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def interval_union(intervals):
    result = []
    for start, end in sorted(intervals):
        if not result or start > result[-1][1]:
            result.append([start, end])
        elif end > result[-1][1]:
            result[-1][1] = end
    return result


class CoverageIndex:
    def __init__(self, features, strand_specific):
        groups = defaultdict(list)
        for seqid, start, end, strand, _ in features:
            key = (seqid, strand) if strand_specific else seqid
            groups[key].append((start, end))
        self.blocks = {k: interval_union(v) for k, v in groups.items()}
        self.ends = {k: [end for _, end in v] for k, v in self.blocks.items()}
        self.strand_specific = strand_specific

    def fraction(self, feature):
        seqid, start, end, strand, _ = feature
        key = (seqid, strand) if self.strand_specific else seqid
        blocks = self.blocks.get(key, [])
        ends = self.ends.get(key, [])
        i = bisect_right(ends, start)
        covered = 0
        while i < len(blocks) and blocks[i][0] < end:
            left, right = blocks[i]
            covered += max(0, min(end, right) - max(start, left))
            i += 1
        return covered / (end - start)


def coverage_counts(queries, targets, strand_specific):
    index = CoverageIndex(targets, strand_specific)
    full = partial = none = 0
    total_fraction = 0.0
    for feature in queries:
        frac = index.fraction(feature)
        total_fraction += frac
        if frac >= 1 - 1e-12:
            full += 1
        elif frac > 0:
            partial += 1
        else:
            none += 1
    return full, partial, none, round(100 * total_fraction / len(queries), 2) if queries else None


def load_gff(path):
    features = defaultdict(list)
    counts = Counter()
    seqids = set()
    meta = {}
    seq_regions = {}
    invalid_rows = 0
    out_of_bounds = 0
    pseudo_cds = 0
    pseudo_genes = 0
    protein_genes = []
    all_feature_ends = []

    with opengz(path) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        seq_regions[parts[1]] = (int(parts[2]), int(parts[3]))
                    except ValueError:
                        pass
            for directive, label in (
                ("#!processor", "processor"),
                ("#!genome-build-accession", "assembly_accession"),
                ("#!genome-build", "assembly_build"),
                ("#!annotation-date", "annotation_date"),
            ):
                if line.startswith(directive + " "):
                    meta[label] = line.strip()[len(directive) + 1:]
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                invalid_rows += 1
                continue
            try:
                start1, end1 = int(f[3]), int(f[4])
            except ValueError:
                invalid_rows += 1
                continue
            if start1 < 1 or end1 < start1:
                invalid_rows += 1
                continue
            seqid, feature_type, strand = f[0], f[2], f[6]
            if seqid in seq_regions:
                lower, upper = seq_regions[seqid]
                if start1 < lower or end1 > upper:
                    out_of_bounds += 1
            attrs = attr_dict(f[8])
            record = (seqid, start1 - 1, end1, strand, attrs)
            seqids.add(seqid)
            counts[feature_type] += 1
            all_feature_ends.append((seqid, end1))
            if feature_type in FEATURES:
                features[feature_type].append(record)
            if feature_type == "gene":
                if attrs.get("gene_biotype", attrs.get("biotype")) == "protein_coding":
                    protein_genes.append(record)
                if attrs.get("pseudo", "").lower() in ("true", "1") or "pseudogene" in attrs:
                    pseudo_genes += 1
            if feature_type == "CDS":
                if attrs.get("pseudo", "").lower() in ("true", "1") or "pseudogene" in attrs:
                    pseudo_cds += 1

    return (features, counts, seqids, meta, invalid_rows, out_of_bounds,
            protein_genes, pseudo_cds, pseudo_genes, all_feature_ends)


def fasta_sizes(path):
    result = {}
    seqid = None
    size = 0
    with opengz(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if seqid is not None:
                    result[seqid] = size
                seqid = line[1:].split()[0]
                size = 0
            else:
                size += len(line.strip())
    if seqid is not None:
        result[seqid] = size
    return result


def audit(folder, check_fasta):
    genome = folder.name
    row = {"genome": genome, "status": "OK"}
    gff = get_file(folder, "genomic.gff")
    if gff is None:
        row["status"] = "MISSING_GFF"
        return row

    (features, counts, seqids, meta, invalid_rows, out_of_bounds,
     protein_genes, pseudo_cds, pseudo_genes, all_feature_ends) = load_gff(gff)
    row.update({"gff": str(gff), "gff_seqids": len(seqids),
                "processor": meta.get("processor", ""),
                "assembly_build": meta.get("assembly_build", ""),
                "assembly_accession": meta.get("assembly_accession", ""),
                "annotation_date": meta.get("annotation_date", ""),
                "invalid_feature_rows": invalid_rows,
                "gff_bounds_violations": out_of_bounds,
                "pseudo_flagged_CDS": pseudo_cds,
                "pseudo_flagged_gene": pseudo_genes,
                "protein_coding_genes": len(protein_genes)})
    for feature in FEATURES:
        row[feature] = counts[feature]

    gene = features["gene"]
    cds = features["CDS"]
    pseudo = features["pseudogene"]
    rna = features["tRNA"] + features["rRNA"] + features["tmRNA"]
    checks = {
        "gene_by_CDS_same": (gene, cds, True),
        "protein_gene_by_CDS_same": (protein_genes, cds, True),
        "gene_by_CDS_or_RNA_same": (gene, cds + rna, True),
        "RNA_by_gene_same": (rna, gene, True),
        "exon_by_RNA_same": (features["exon"], rna, True),
        "CDS_by_gene_same": (cds, gene, True),
        "CDS_by_gene_any": (cds, gene, False),
        "pseudo_by_CDS_same": (pseudo, cds, True),
        "pseudo_by_CDS_any": (pseudo, cds, False),
        "pseudo_by_gene_any": (pseudo, gene, False),
    }
    for name, (query, target, same_strand) in checks.items():
        full, partial, none, mean = coverage_counts(query, target, same_strand)
        row[f"{name}_full"] = full
        row[f"{name}_partial"] = partial
        row[f"{name}_none"] = none
        row[f"{name}_mean_pct"] = "" if mean is None else mean

    if check_fasta:
        fasta = get_file(folder, "genomic.fna")
        if fasta is None:
            row["fasta_status"] = "MISSING_FASTA"
        else:
            sizes = fasta_sizes(fasta)
            row["fasta_status"] = "OK"
            row["gff_seqids_missing_from_fasta"] = len(seqids - sizes.keys())
            row["fasta_seqids_missing_from_gff"] = len(sizes.keys() - seqids)
            row["features_beyond_fasta_length"] = sum(
                end > sizes[seqid]
                for seqid, end in all_feature_ends
                if seqid in sizes
            )
    if invalid_rows or out_of_bounds:
        row["status"] = "INVALID_GFF_ROWS"
    elif check_fasta and row.get("fasta_status") != "OK":
        row["status"] = "MISSING_FASTA"
    elif check_fasta and (row.get("gff_seqids_missing_from_fasta", 0) > 0 or
                          row.get("fasta_seqids_missing_from_gff", 0) > 0 or
                          row.get("features_beyond_fasta_length", 0) > 0):
        row["status"] = "GFF_FASTA_MISMATCH"
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, required=True,
                    help="Parent directory containing GCA_*/genomic.gff")
    ap.add_argument("--output", type=Path, default=Path("annotation_audit_162.tsv"))
    ap.add_argument("--check-fasta", action="store_true", help="Read genomic.fna to check seqIDs and lengths")
    args = ap.parse_args()
    if not args.root.is_dir():
        ap.error(f"Root directory does not exist: {args.root}")
    folders = sorted(p for p in args.root.iterdir() if p.is_dir() and p.name.startswith(("GCA_", "GCF_")))
    if not folders:
        ap.error(f"No GCA_ or GCF_ directories found under: {args.root}")
    rows = []
    for i, folder in enumerate(folders, 1):
        try:
            row = audit(folder, args.check_fasta)
        except Exception as exc:
            row = {"genome": folder.name, "status": f"ERROR: {type(exc).__name__}: {exc}"}
        rows.append(row)
        print(f"[{i:>3}/{len(folders)}] {folder.name}: {row['status']}", flush=True)

    headers = ["genome", "status"] + sorted(k for row in rows for k in row if k not in ("genome", "status"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    ok = sum(r["status"] == "OK" for r in rows)
    print(f"\nAudited {len(rows)} assembly directories ({ok} structurally OK).")
    print(f"Full report: {args.output}")
    with_pseudo = [r for r in rows if int(r.get("pseudogene", 0)) > 0]
    print(f"GFFs with explicit 'pseudogene' features: {len(with_pseudo)}/{len(rows)}")
    without_cds = [r for r in with_pseudo if int(r.get("pseudo_by_CDS_any_none", 0)) > 0]
    print(f"GFFs with pseudogene loci not covered by any CDS: {len(without_cds)}")
    print("Note: missing gene/CDS pairings or zero pseudogene labels do not by themselves prove incomplete biological annotation.")

if __name__ == "__main__":
    main()
