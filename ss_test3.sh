#!/bin/bash

# === CONFIGURATION ===
# This script follows the instructions by re-aligning the pair of sequence
# fragments found in the tblastn output using the ssearch36 tool.

#------------SSEARCH-------------------
echo "Checking and downloading SSEARCH..."
if [ ! -f fasta36-linux64.tar.gz ]; then
    wget https://fasta.bioch.virginia.edu/wrpearson/fasta/CURRENT/fasta36-linux64.tar.gz 
    tar -xzf fasta36-linux64.tar.gz 
else
    echo "SSEARCH archive already exists. Skipping download."
fi

echo "--- Filtering results based on E-value threshold ---"

# Check if the filtered output file does NOT already exist (-f checks for a file)
if [ ! -f "genome.tblastn.filtered.out" ]; then
    echo "Filtered file not found. Running the filtering process..."
    
    # Run the time-consuming awk command
    awk -v threshold="$MIN_EVALUE" '$11 <= threshold' genome.tblastn.out > genome.tblastn.filtered.out
    
    echo "Filtering complete. Results saved to genome.tblastn.filtered.out"
else
    # If the file exists, print a message and skip the step
    echo "Filtered file 'genome.tblastn.filtered.out' already exists. Skipping."
fi

echo "Filtering complete. Results saved to genome.tblastn.filtered.out"
# You can see how many significant hits were found with:
wc -l genome.tblastn.filtered.out

echo "Original size"
wc -l genome.tblastn.out

# Automatically determine the number of available CPU cores
total_cores=$(nproc)

# Set NUM_JOBS to be the total number of cores minus one.
if (( total_cores > 1 )); then
    NUM_JOBS=$((total_cores - 1))
else
    NUM_JOBS=1
fi

# Input and output files
TBLASTN_OUTPUT="genome.tblastn.filtered.out"
FINAL_RESULTS="ssearch_fragment_realignment_results.out"
SSEARCH_BIN="fasta-36.3.8i/bin/ssearch36"
# =====================


# --- Pre-flight checks ---
echo "Checking for necessary files..."
for f in "$TBLASTN_OUTPUT" "$SSEARCH_BIN"; do
    if [ ! -f "$f" ]; then
        echo "Error: Required file not found: $f"
        exit 1
    fi
done

# Clean previous results file
> "$FINAL_RESULTS"

# --- Define the core processing function ---
# This version aligns the two sequence fragments from the tblastn output.
process_fragments() {
    local qseq_fragment=$1
    local sseq_fragment=$2
    local ssearch_bin=$3
    local line_num=$4 # Use line number for reference

    # Clean the fragments by removing gap characters ('-')
    local clean_qseq=${qseq_fragment//-/}
    local clean_sseq=${sseq_fragment//-/}

    # Execute ssearch using process substitution to avoid temp files
    local alignment_output
    alignment_output=$("$ssearch_bin" -b 1 -d 1 -k 500 -f -11 -g -1 -m 10 -q -s BL62 -z 11 -Z 10000 \
        <(echo -e ">query_frag_${line_num}\n${clean_qseq}") \
        <(echo -e ">subject_frag_${line_num}\n${clean_sseq}")
    )
    
    # Use the validated filter to check for success
    if grep -q "The best scores are:" <<< "$alignment_output"; then
        # Add the original line number to the output for reference
        echo "# Original TBLASTN hit line: $line_num"
        echo "$alignment_output"
    fi
}

export -f process_fragments

# --- Main Execution using GNU Parallel ---
echo "Starting parallel re-alignment of TBLASTN fragments with $NUM_JOBS jobs..."
echo "This will run on the filtered dataset."

# Use awk to extract the two fragment columns (9 and 10)
awk '{print $9, $10}' "$TBLASTN_OUTPUT" | \
    parallel --eta -j "$NUM_JOBS" --colsep ' ' --no-notice --line-buffer --link \
    "process_fragments {1} {2} $SSEARCH_BIN {#}" >> "$FINAL_RESULTS"


echo "---"
echo "Parallel processing complete."
echo "All significant fragment alignment results have been combined into $FINAL_RESULTS"
