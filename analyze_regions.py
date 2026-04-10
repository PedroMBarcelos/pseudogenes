#!/usr/bin/env python3

import sys
import re
from collections import defaultdict
import argparse

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


def find_disablements(query_seq, subject_seq):
    """
    Scans a pairwise alignment for frameshifts and internal stop codons.

    Args:
        query_seq (str): The reference protein sequence from the alignment.
        subject_seq (str): The translated genomic sequence from the alignment.

    Returns:
        tuple: (frameshift_count, stop_codon_count)
    """
    frameshifts = 0
    internal_stops = 0

    # Count internal stop codons ('*') in the subject sequence
    internal_stops = subject_seq[:-1].count('*')

    # Find frameshifts by looking at gaps ('-')
    gap_pattern = re.compile(r"(-+)")
    
    for gap_match in gap_pattern.finditer(subject_seq):
        if len(gap_match.group(1)) % 3 != 0:
            frameshifts += 1
            
    for gap_match in gap_pattern.finditer(query_seq):
        if len(gap_match.group(1)) % 3 != 0:
            frameshifts += 1

    return frameshifts, internal_stops


def main(args):
    """Main execution function."""
    
    ssearch_data = parse_ssearch_results(args.ssearch, args.min_evalue)

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
            # Add a more robust check for the merged regions file format.
            if len(parts) < 4:
                print(f"FATAL: Malformed line in merged regions file '{args.merged_regions}'.", file=sys.stderr)
                print(f"       Expected 4 columns (chrom, start, end, strand), but found {len(parts)}.", file=sys.stderr)
                print(f"       Line content: {line}", file=sys.stderr)
                print("HINT:  Ensure the bedtools merge command in your shell script includes '-c 6 -o first' to output the strand.", file=sys.stderr)
                sys.exit(1)
            
            region_chrom, region_start, region_end, region_strand = parts[0], int(parts[1]), int(parts[2]), parts[3]
            region_id = f"region_{i+1}"

            proteins_in_region = defaultdict(lambda: {'total_score': 0, 'total_length': 0, 'hits': []})

            for hit in all_hits.get(region_chrom, []):
                if hit['strand'] == region_strand and hit['start'] < region_end and hit['end'] > region_start:
                    line_num = hit['line_num']
                    if line_num in ssearch_data:
                        qseqid = hit['qseqid']
                        ssearch_hit = ssearch_data[line_num]
                        proteins_in_region[qseqid]['total_score'] += ssearch_hit.get('score', 0)
                        proteins_in_region[qseqid]['total_length'] += ssearch_hit.get('length', 0)
                        proteins_in_region[qseqid]['hits'].append(ssearch_hit)
            
            if not proteins_in_region:
                continue

            best_protein = None
            max_density = -1.0

            for qseqid, data in proteins_in_region.items():
                if data['total_length'] > 0:
                    density = data['total_score'] / data['total_length']
                    if density > max_density:
                        max_density = density
                        best_protein = qseqid
            
            if not best_protein:
                continue

            total_frameshifts = 0
            total_stops = 0
            parent_alignments = proteins_in_region[best_protein]['hits']
            
            for alignment in parent_alignments:
                frameshifts, stops = find_disablements(alignment['q_seq'], alignment['s_seq'])
                total_frameshifts += frameshifts
                total_stops += stops
            
            total_disablements = total_frameshifts + total_stops

            final_results.append({
                'RegionID': region_id,
                'Chromosome': region_chrom,
                'Start': region_start,
                'End': region_end,
                'Strand': region_strand,
                'ParentProtein': best_protein,
                'ScoreDensity': f"{max_density:.4f}",
                'NumHSPs': len(parent_alignments),
                'TotalFrameshifts': total_frameshifts,
                'InternalStops': total_stops,
                'TotalDisablements': total_disablements
            })

    print(f"INFO: Writing final report to {args.output}...")
    with open(args.output, 'w') as f_out:
        f_out.write("\t".join(report_columns) + "\n")
        for row in final_results:
            f_out.write("\t".join(str(row[column]) for column in report_columns) + "\n")

    if not final_results:
        print("WARNING: No results were generated.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyze merged genomic regions to find parent proteins and pseudogene features.")
    parser.add_argument("--ssearch", required=True, help="Path to the concatenated ssearch36 results file.")
    parser.add_argument("--bed_hits", required=True, help="Path to the sorted BED file of all original TBLASTN HSPs.")
    parser.add_argument("--merged_regions", required=True, help="Path to the BED file of merged genomic regions.")
    parser.add_argument("--output", required=True, help="Path for the final TSV report.")
    parser.add_argument("--min_evalue", required=True, type=float, help="E-value threshold for filtering SSEARCH results.")
    
    args = parser.parse_args()
    main(args)
