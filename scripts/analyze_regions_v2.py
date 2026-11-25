#!/usr/bin/env python3

import sys
import re
from collections import defaultdict
import argparse
import csv

from intervaltree import Interval, IntervalTree
from utils import select_best_protein

"""
This script performs an advanced analysis of genomic regions identified from TBLASTN and SSEARCH results.

Workflow:
1.  Parse SSEARCH alignment results in the '-m 10' machine-readable format, filtering them by a significant E-value threshold.
2.  For each merged genomic region (from BEDTools), identify all the original HSPs that fall within it.
3.  For each UniProt protein hitting a merged region, calculate its "score density" as the sum of its HSP scores divided by the sum of its HSP lengths.
4.  Select the UniProt protein with the highest score density as the "parent protein" for that region.
5.  Re-examine the alignments corresponding to the parent protein to find potential disablements (frameshifts, internal stop codons).
6.  Generate a final report summarizing the findings for each merged region.
"""

def parse_tblastn_results(tblastn_file):
    tblastn_data = {}

    with open(tblastn_file, 'r') as f:
        for line_num, line in enumerate(f, start=1):
            if not line.strip():
                continue

            fields = line.rstrip('\n').split('\t')

            if len(fields) < 25:
                print(
                    f"WARNING: Skipping malformed TBLASTN "
                    f"line {line_num}: expected 25 columns, "
                    f"found {len(fields)}",
                    file=sys.stderr
                )
                continue

            tblastn_data[line_num] = {
                "qseqid": fields[0],
                "sseqid": fields[2],
                "qseq": fields[8],
                "sseq": fields[9],
                "evalue": float(fields[10]),
                "sstart": int(fields[6]),
                "send": int(fields[7]),
                "sframe": int(fields[21]),
                "sstrand": fields[22]
            }

    print(
        f"INFO: Parsed {len(tblastn_data)} "
        f"TBLASTN alignments.",
        file=sys.stderr
    )

    return tblastn_data


def parse_ssearch_results(ssearch_file, min_evalue_threshold):
    """
    Parses the concatenated ssearch36 output file in the '-m 10' format,
    applying an E-value filter.

    Args:
        ssearch_file (str): Path to the ssearch results file.
        min_evalue_threshold (float): The E-value cutoff. SSEARCH hits with an
                                      E-value greater than this will be ignored.

    Returns:
        dict: A dictionary mapping the original TBLASTN line number to a dict
              containing score, length, and alignment details.
    """
    print(f"INFO: Parsing SSEARCH '-m 10' results with E-value threshold <= {min_evalue_threshold}...")
    ssearch_data = {}
    current_entry = {}
    capture_state = None  # Can be 'query', 'subject', or None

    with open(ssearch_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # A new record starts with our custom comment
            if line.startswith("# Original TBLASTN hit line:"):
                # If we have a completed entry, process and store it
                if current_entry and 'line_num' in current_entry:
                    # Finalize sequence assembly
                    current_entry['q_seq'] = "".join(current_entry.get('q_seq_parts', []))
                    current_entry['s_seq'] = "".join(current_entry.get('s_seq_parts', []))
                    
                    # Filter and store
                    if current_entry.get('evalue', float('inf')) <= min_evalue_threshold:
                        ssearch_data[current_entry['line_num']] = current_entry

                # Start a new entry
                current_entry = {
                    'line_num': int(line.split()[-1]),
                    'q_seq_parts': [],
                    's_seq_parts': []
                }
                capture_state = None
                continue

            if not current_entry:
                continue

            # Extract key-value data
            if line.startswith('; sw_expect:'):
                current_entry['evalue'] = float(line.split()[-1])
            elif line.startswith('; sw_score:'):
                current_entry['score'] = int(line.split()[-1])
            elif line.startswith('; sq_len:') and capture_state == 'query':
                current_entry['length'] = int(line.split()[-1])
            
            # Control sequence capture state
            elif line.startswith('>query_frag_'):
                capture_state = 'query'
            elif line.startswith('>>/dev/fd/62'):
                capture_state = 'subject'
            elif line.startswith('; al_cons:'):
                capture_state = None # End of alignment block
            
            # Capture sequence data based on state
            elif capture_state == 'query' and not line.startswith((';', '>')):
                current_entry['q_seq_parts'].append(line)
            elif capture_state == 'subject' and not line.startswith((';', '>>')):
                current_entry['s_seq_parts'].append(line)

    # Process the very last entry in the file
    if current_entry and 'line_num' in current_entry:
        current_entry['q_seq'] = "".join(current_entry.get('q_seq_parts', []))
        current_entry['s_seq'] = "".join(current_entry.get('s_seq_parts', []))
        if current_entry.get('evalue', float('inf')) <= min_evalue_threshold:
            ssearch_data[current_entry['line_num']] = current_entry

    print(f"INFO: Parsed and kept {len(ssearch_data)} significant SSEARCH alignments.")
    return ssearch_data


def count_internal_stops(subject_seq):
    """
    Count stop symbols explicitly present in the original TBLASTN
    translated subject alignment.

    This does not infer genomic frameshifts. A stop is counted only
    if TBLASTN actually includes '*' in its sseq field.
    """
    cleaned_subject = subject_seq.rstrip("-")

    if not cleaned_subject:
        return 0

    total_stops = cleaned_subject.count("*")

    # Preserve the existing convention of excluding a terminal stop.
    if cleaned_subject.endswith("*"):
        total_stops -= 1

    return max(0, total_stops)

def main(args):
    """Main execution function."""
    
    ssearch_data = parse_ssearch_results(args.ssearch, args.min_evalue)
    tblastn_data = parse_tblastn_results(args.tblastn)

    report_columns = [
        'RegionID',
        'Chromosome',
        'Start',
        'End',
        'Strand',
        'ParentProtein',
        'ScoreDensity',
        'NumHSPs',
        'TotalFrameshifts',
        'InternalStops',
        'TotalDisablements'
    ]

    print("INFO: Loading all HSP hits from BED file...")
    all_hits = defaultdict(list)
    with open(args.bed_hits, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7: continue
            chrom, start, end, qseqid, score, strand, line_num = parts
            all_hits[chrom].append({
                'start': int(start), 'end': int(end), 'qseqid': qseqid,
                'score': float(score), 'strand': strand, 'line_num': int(line_num)
            })

    print("INFO: Analyzing each merged region...")
    final_results = []

    with open(args.merged_regions, 'r') as f_merged:
        for i, line in enumerate(f_merged):
            line = line.strip()
            if not line: continue

            parts = line.split('\t')

            if len(parts) < 4:
                print(
                    f"FATAL: Malformed line in merged regions file "
                    f"'{args.merged_regions}'.",
                    file=sys.stderr
                )
                # FIX 1: Indent these lines so they are inside the `if` block
                print(
                    f"Expected 4 columns, found {len(parts)}.",
                    file=sys.stderr
                )
                sys.exit(1)

            region_chrom = parts[0]
            region_start = int(parts[1])
            region_end = int(parts[2])
            region_strand = parts[3]
            region_id = f"region_{i+1}"

            # FIX 2: Indent everything from here down to final_results.append 
            # by an additional 4 spaces so it remains inside the `for` loop.
            proteins_in_region = defaultdict(
                lambda: {
                    'total_score': 0,
                    'total_length': 0,
                    'hits': []
                }
            )

            for hit in all_hits.get(region_chrom, []):
                if (
                    hit['strand'] == region_strand
                    and hit['start'] < region_end
                    and hit['end'] > region_start
                ):
                    line_num = hit['line_num']

                    if line_num in ssearch_data:
                        qseqid = hit['qseqid']
                        ssearch_hit = ssearch_data[line_num]

                        proteins_in_region[qseqid]['total_score'] += (
                            ssearch_hit.get('score', 0)
                        )

                        proteins_in_region[qseqid]['total_length'] += (
                            ssearch_hit.get('length', 0)
                        )

                        # NEW: Preserve the original TBLASTN row number.
                        proteins_in_region[qseqid]['hits'].append({
                            **ssearch_hit,
                            'line_num': line_num
                        })

            if not proteins_in_region:
                continue

            best_protein = select_best_protein(proteins_in_region)
            if not best_protein:
                continue

            p_data = proteins_in_region[best_protein]

            total_score = p_data['total_score']
            total_length = p_data['total_length']

            score_density = (
                total_score / total_length
                if total_length > 0 else 0.0
            )

            parent_alignments = p_data['hits']

            total_stops = 0
            seen_tblastn_lines = set()

            for alignment in parent_alignments:
                line_num = alignment['line_num']

                if line_num in seen_tblastn_lines:
                    continue

                seen_tblastn_lines.add(line_num)

                original_hit = tblastn_data.get(line_num)

                if original_hit is None:
                    print(
                        f"WARNING: Missing TBLASTN alignment "
                        f"for line {line_num}",
                        file=sys.stderr
                    )
                    continue

                total_stops += count_internal_stops(
                    original_hit['sseq']
                )

            # Frameshifts cannot be determined reliably from
            # the current SSEARCH protein alignments.
            total_frameshifts = 0
            total_disablements = total_stops + total_frameshifts

            final_results.append({
                'RegionID': region_id,
                'Chromosome': region_chrom,
                'Start': region_start,
                'End': region_end,
                'Strand': region_strand,
                'ParentProtein': best_protein,
                'ScoreDensity': f"{score_density:.4f}",
                'NumHSPs': len(parent_alignments),
                'TotalFrameshifts': total_frameshifts,
                'InternalStops': total_stops,
                'TotalDisablements': total_disablements
            })

    print(f"INFO: Writing final report to {args.output}...")
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=report_columns,
            delimiter='\t'
        )

        writer.writeheader()

        if final_results:
            writer.writerows(final_results)

    print(
        f"INFO: Successfully wrote {len(final_results)} regions "
        f"to {args.output}"
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyze merged genomic regions to find parent proteins and pseudogene features.")
    parser.add_argument("--ssearch", required=True, help="Path to the concatenated ssearch36 results file.")
    parser.add_argument("--bed_hits", required=True, help="Path to the sorted BED file of all original TBLASTN HSPs.")
    parser.add_argument("--merged_regions", required=True, help="Path to the BED file of merged genomic regions.")
    parser.add_argument("--output", required=True, help="Path for the final TSV report.")
    parser.add_argument("--min_evalue", required=True, type=float, help="E-value threshold for filtering SSEARCH results.")
    parser.add_argument("--tblastn", required=True, help="Path to the filtered original TBLASTN tabular output.")
    
    args = parser.parse_args()
    main(args)
