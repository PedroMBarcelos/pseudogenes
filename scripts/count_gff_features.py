#!/usr/bin/env python3
"""Count GFF3 feature TYPE labels across an assembly collection (not loci).

Processes genomic.gff or genomic.gff.gz once per GCA_/GCF_ directory. Stops at
##FASTA and ignores comment/directive rows. Counts 'gene' and 'CDS' as distinct
features, even when they describe the same locus.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--totals", type=Path, default=Path("gff_feature_counts_total.tsv"))
    ap.add_argument("--by-genome", type=Path, default=Path("gff_feature_counts_by_genome.tsv"))
    args = ap.parse_args()
    if not args.root.is_dir():
        ap.error(f"Directory not found: {args.root}")

    overall = Counter()
    each = []
    missing = []
    bad = []
    folders = sorted(p for p in args.root.iterdir()
                     if p.is_dir() and p.name.startswith(("GCA_", "GCF_")))
    for folder in folders:
        gff = next((p for p in (folder / "genomic.gff", folder / "genomic.gff.gz")
                    if p.is_file()), None)
        if gff is None:
            missing.append(folder.name)
            continue
        counts = Counter()
        with open_text(gff) as handle:
            for line_number, line in enumerate(handle, 1):
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    bad.append((folder.name, line_number))
                    continue
                counts[fields[2]] += 1
        overall.update(counts)
        each.append((folder.name, counts))

    args.totals.parent.mkdir(parents=True, exist_ok=True)
    with args.totals.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(("feature_type", "total_features", "genomes_with_feature"))
        for ftype, n in overall.most_common():
            writer.writerow((ftype, n, sum(ftype in c for _, c in each)))

    args.by_genome.parent.mkdir(parents=True, exist_ok=True)
    with args.by_genome.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(("genome", "feature_type", "count"))
        for genome, counts in each:
            for feature, n in counts.most_common():
                writer.writerow((genome, feature, n))

    print(f"Assemblies: {len(folders)}; GFFs: {len(each)}; missing GFF: {len(missing)}")
    print(f"Malformed GFF rows: {len(bad)}")
    print(f"\n{'Feature':24} {'Total':>10} {'Genomes':>10}")
    for ftype, n in overall.most_common():
        print(f"{ftype:24} {n:10,d} {sum(ftype in c for _, c in each):10,d}")
    print(f"\nTotals: {args.totals}\nPer genome: {args.by_genome}")


if __name__ == "__main__":
    main()
