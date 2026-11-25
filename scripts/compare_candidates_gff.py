#!/usr/bin/env python3
"""Compare protein-guided candidate regions (BED coordinates in TSV) with GFF3 features.

Exclude GFF3 feature type 'region'. Preserve all candidate TSV rows in the
summary and generate one detail row per candidate–GFF feature overlap.
"""
import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

INCLUDED_FEATURE_TYPES = {
    "gene",
    "CDS",
    "pseudogene",
    "exon",
    "tRNA",
    "rRNA",
    "tmRNA",
}

def open_text(path):
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8') if path.suffix == '.gz' else path.open(encoding='utf-8')


def parse_attributes(field):
    result = {}
    for entry in field.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            result[key] = unquote(value)
    return result


def load_gff(path):
    by_chrom = defaultdict(list)
    count = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{line_number}: expected 9 GFF columns, got {len(fields)}')
            seqid, source, feature_type, begin, finish, score, strand, phase, attr_text = fields
            if feature_type not in INCLUDED_FEATURE_TYPES:
                continue
            start_1, end_1 = int(begin), int(finish)
            if start_1 < 1 or end_1 < start_1:
                raise ValueError(f'{path}:{line_number}: invalid GFF coordinates')
            attrs = parse_attributes(attr_text)
            by_chrom[seqid].append({
                'start0': start_1 - 1,
                'end0': end_1,
                'start1': start_1,
                'end1': end_1,
                'type': feature_type,
                'strand': strand,
                'source': source,
                'id': attrs.get('ID', ''),
                'parent': attrs.get('Parent', ''),
                'name': attrs.get('Name', attrs.get('gene', '')),
                'gene': attrs.get('gene', ''),
                'locus_tag': attrs.get('locus_tag', ''),
                'pseudo': attrs.get('pseudo', ''),
                'product': attrs.get('product', ''),
                'note': attrs.get('Note', ''),
                'line': line_number,
            })
            count += 1
    for features in by_chrom.values():
        features.sort(key=lambda x: (x['start0'], x['end0']))
    return by_chrom, count


def merged_coverage(intervals):
    """Union of intersected BED intervals, without double-counting overlapping GFF features."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    begin, end = intervals[0]
    total = 0
    for a, b in intervals[1:]:
        if a <= end:
            end = max(end, b)
        else:
            total += end - begin
            begin, end = a, b
    return total + end - begin


def find_overlaps(row, by_chrom):
    chromosome = row['chromosome']
    start = int(row['start'])
    end = int(row['end'])
    if start < 0 or end <= start:
        raise ValueError(f"Invalid candidate BED interval: {row.get('region_id')}: {start}-{end}")
    for feature in by_chrom.get(chromosome, []):
        if feature['start0'] >= end:
            break
        left = max(start, feature['start0'])
        right = min(end, feature['end0'])
        if left < right:
            yield feature, left, right


def strand_relation(candidate_strand, feature_strand):
    if candidate_strand not in ('+', '-') or feature_strand not in ('+', '-'):
        return 'unknown'
    return 'same' if candidate_strand == feature_strand else 'opposite'


def compare(candidate_tsv, gff_path, outdir):
    # Long genomic DNA / JSON cells may exceed the csv module's default field limit.
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            limit //= 10

    by_chrom, n_gff = load_gff(gff_path)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = candidate_tsv.stem
    summary_path = outdir / f'{prefix}_gff_summary.tsv'
    details_path = outdir / f'{prefix}_gff_overlaps.tsv'

    extra_columns = [
        'gff_overlap_status', 'gff_n_features', 'gff_feature_types',
        'gff_feature_ids', 'gff_gene_names', 'gff_strand_relations',
        'gff_covered_nt_union', 'gff_covered_candidate_pct',
    ]
    details_columns = [
        'region_id', 'chromosome', 'candidate_start0', 'candidate_end0',
        'candidate_strand', 'classification', 'protein_id',
        'gff_type', 'gff_id', 'gff_name', 'gff_parent', 'gff_gene',
        'gff_locus_tag', 'gff_pseudo', 'gff_product', 'gff_note', 'gff_source', 'gff_start1', 'gff_end1',
        'gff_strand', 'strand_relation', 'overlap_nt',
        'candidate_overlap_pct', 'gff_feature_overlap_pct',
    ]

    n_candidates = 0
    n_overlapping = 0
    n_pairs = 0
    missing_seqids = set()
    feature_type_counts = Counter()
    with candidate_tsv.open(encoding='utf-8', newline='') as inp, \
         summary_path.open('w', encoding='utf-8', newline='') as summary_handle, \
         details_path.open('w', encoding='utf-8', newline='') as details_handle:
        reader = csv.DictReader(inp, delimiter='\t')
        needed = {'region_id', 'chromosome', 'start', 'end', 'strand'}
        if not reader.fieldnames or not needed.issubset(reader.fieldnames):
            raise ValueError(f'Candidate TSV missing columns: {needed - set(reader.fieldnames or [])}')
        summary_writer = csv.DictWriter(summary_handle, fieldnames=reader.fieldnames + extra_columns, delimiter='\t')
        details_writer = csv.DictWriter(details_handle, fieldnames=details_columns, delimiter='\t')
        summary_writer.writeheader()
        details_writer.writeheader()

        for row in reader:
            n_candidates += 1
            chrom = row['chromosome']
            start, end = int(row['start']), int(row['end'])
            if chrom not in by_chrom:
                missing_seqids.add(chrom)
            hits = list(find_overlaps(row, by_chrom))
            if hits:
                n_overlapping += 1
            type_counter = Counter()
            feature_ids = set()
            gene_names = set()
            strand_relations = set()
            intervals = []

            for feature, left, right in hits:
                overlap = right - left
                relation = strand_relation(row['strand'], feature['strand'])
                type_counter[feature['type']] += 1
                feature_type_counts[feature['type']] += 1
                n_pairs += 1
                intervals.append((left, right))
                strand_relations.add(relation)
                if feature['id']:
                    feature_ids.add(feature['id'])
                if feature['type'] in {'gene', 'pseudogene'}:
                    label = feature['gene'] or feature['locus_tag'] or feature['name'] or feature['id']
                    if label:
                        gene_names.add(label)
                details_writer.writerow({
                    'region_id': row['region_id'], 'chromosome': chrom,
                    'candidate_start0': start, 'candidate_end0': end,
                    'candidate_strand': row['strand'],
                    'classification': row.get('region_type', row.get('classification', '')),
                    'protein_id': row.get('protein_id', ''),
                    'gff_type': feature['type'], 'gff_id': feature['id'],
                    'gff_name': feature['name'], 'gff_parent': feature['parent'],
                    'gff_gene': feature['gene'], 'gff_locus_tag': feature['locus_tag'],
                    'gff_pseudo': feature['pseudo'], 'gff_product': feature['product'],
                    'gff_note': feature['note'], 'gff_source': feature['source'],
                    'gff_start1': feature['start1'], 'gff_end1': feature['end1'],
                    'gff_strand': feature['strand'], 'strand_relation': relation,
                    'overlap_nt': overlap,
                    'candidate_overlap_pct': round(100 * overlap / (end - start), 3),
                    'gff_feature_overlap_pct': round(100 * overlap / (feature['end0'] - feature['start0']), 3),
                })

            covered = merged_coverage(intervals)
            row.update({
                'gff_overlap_status': 'overlap' if hits else 'no_overlap',
                'gff_n_features': len(hits),
                'gff_feature_types': ';'.join(f'{t}:{n}' for t, n in sorted(type_counter.items())),
                'gff_feature_ids': ';'.join(sorted(feature_ids)),
                'gff_gene_names': ';'.join(sorted(gene_names)),
                'gff_strand_relations': ';'.join(sorted(strand_relations)),
                'gff_covered_nt_union': covered,
                'gff_covered_candidate_pct': round(100 * covered / (end - start), 3),
            })
            summary_writer.writerow(row)

    print(f'GFF features indexed (excluding type region): {n_gff}')
    print(f'Candidates: {n_candidates}')
    print(f'Candidates overlapping any GFF feature: {n_overlapping}')
    print(f'Candidates without GFF overlap: {n_candidates - n_overlapping}')
    print(f'Candidate–feature overlap pairs: {n_pairs}')
    print(f'Feature types in overlap pairs: {dict(sorted(feature_type_counts.items()))}')
    if missing_seqids:
        print('WARNING: candidate chromosome IDs absent from GFF:', ', '.join(sorted(missing_seqids)), file=sys.stderr)
    print(f'Summary: {summary_path}')
    print(f'Per-feature overlaps: {details_path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', type=Path, required=True, help='Final TSV, with BED start/end coordinates')
    parser.add_argument('--gff', type=Path, required=True, help='Genome GFF3 (or .gz)')
    parser.add_argument('--outdir', type=Path, required=True, help='Output directory')
    args = parser.parse_args()
    compare(args.candidates, args.gff, args.outdir)


if __name__ == '__main__':
    main()
