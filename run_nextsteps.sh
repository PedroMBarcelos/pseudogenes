#!/bin/bash

# ==============================================================================
# MASTER ANALYSIS SCRIPT
# ==============================================================================
# This script orchestrates the pipeline:
# 1. Converts TBLASTN output to BED format.
# 2. Sorts the BED file.
# 3. Merges overlapping hits using BEDTools to define genomic regions.
# 4. Calls a Python script to perform the advanced analysis:
#    - Calculate score densities.
#    - Identify a parent protein for each region.
#    - Scan for pseudogene-like disablements.
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status.

# --- Configuration ---
TBLASTN_FILTERED="genome.tblastn.filtered.out"
SSEARCH_RESULTS="ssearch_fragment_realignment_results.out"
PYTHON_SCRIPT="analyze_regions.py"
MIN_EVALUE=$(grep -v '^#' uniprot_sprot_shuffle.blastp.out | sort -k11,11g | head -n 1 | awk '{print $11}')

# --- Output Files ---
BED_HITS="protein_hits.bed"
BED_SORTED="protein_hits.sorted.bed"
BED_MERGED="merged_protein_features.bed"
FINAL_REPORT="pseudogene_analysis_report.tsv"


# --- Pre-flight Checks ---
echo "INFO: Checking for necessary files and programs..."
for f in "$TBLASTN_FILTERED" "$SSEARCH_RESULTS"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required input file not found: $f"
        exit 1
    fi
done

if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools could not be found. Please install it."
    exit 1
fi

if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 could not be found. Please install it."
    exit 1
fi
echo "INFO: All checks passed."


# ==============================================================================
# STEP 1: Convert TBLASTN output to BED format
# ==============================================================================
echo "STEP 1: Converting TBLASTN results to BED format..."
awk 'BEGIN {OFS="\t"} {
    # BED format is 0-based, start < end.
    if ($7 < $8) {
        start = $7 - 1; end = $8; strand = "+";
    } else {
        start = $8 - 1; end = $7; strand = "-";
    }
    # Print: chrom, start, end, name (qseqid), score (bitscore), strand, original_line_number
    # We add the line number to link back to the ssearch results.
    print $3, start, end, $1, $12, strand, NR;
}' "$TBLASTN_FILTERED" > "$BED_HITS"
echo "INFO: Created $BED_HITS"


# ==============================================================================
# STEP 2: Sort the BED file
# ==============================================================================
echo "STEP 2: Sorting BED file for BEDTools..."
sort -k1,1 -k2,2n "$BED_HITS" > "$BED_SORTED"
echo "INFO: Created $BED_SORTED"


# ==============================================================================
# STEP 3: Merge features with BEDTools
# ==============================================================================
echo "STEP 3: Merging features with BEDTools..."
# The `bedtools intersect` command will be used within the Python script
# to associate original hits with the newly defined merged regions.
# Here, we just create the master list of merged regions.
bedtools merge -s -i "$BED_SORTED" -c 6 -o first > "$BED_MERGED"
echo "INFO: Created merged regions file: $BED_MERGED"


# ==============================================================================
# STEP 4: Run analysis with Python
# ==============================================================================
echo "STEP 4: Starting analysis with Python script..."
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: The Python script '$PYTHON_SCRIPT' was not found in this directory."
    exit 1
fi

python3 "$PYTHON_SCRIPT" \
    --ssearch "$SSEARCH_RESULTS" \
    --bed_hits "$BED_SORTED" \
    --merged_regions "$BED_MERGED" \
    --output "$FINAL_REPORT"\
    --min_evalue "$MIN_EVALUE"

echo "=============================================================================="
echo "ANALYSIS COMPLETE!"
echo "The final report on potential pseudogenes is located in: $FINAL_REPORT"
echo "=============================================================================="
