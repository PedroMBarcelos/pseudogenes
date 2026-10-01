import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def analyze_blast_output(file_path="uniprot_sprot_shuffle.blastp.out"):
    """
    Analyzes BLAST output to generate histograms for E-values and bitscores,
    and identifies the minimum E-value and maximum bitscore.
    The script will only proceed if the specified file_path exists.

    Args:
        file_path (str): Path to the BLAST output file (format 6).
    """

    column_names = [
        'qseqid', 'qlen', 'sseqid', 'slen', 'qstart', 'qend', 'sstart', 'send',
        'qseq', 'sseq', 'evalue', 'bitscore', 'score', 'length', 'pident',
        'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos',
        'qcovs', 'qcovhsp'
    ]

    # Check if the BLAST output file exists
    if not os.path.exists(file_path):
        print(f"ERROR: File '{file_path}' not found. Please ensure it exists in the directory.")
        print("Analysis will not proceed without the input file.")
        return # Exit the function if file does not exist

    # Proceed with analysis if the file exists
    print(f"File '{file_path}' found. Proceeding with analysis...")
    try:
        df = pd.read_csv(file_path, sep='\t', header=None, names=column_names)
    except Exception as e: # Catch other potential errors during file reading
        print(f"ERROR reading the file '{file_path}': {e}")
        return

    if df.empty:
        print("The data file is empty. No analysis can be performed.")
        return

    # Convert to numeric, coercing errors for robustness
    df['evalue'] = pd.to_numeric(df['evalue'], errors='coerce')
    df['bitscore'] = pd.to_numeric(df['bitscore'], errors='coerce')

    # Drop rows where conversion failed for key columns
    df.dropna(subset=['evalue', 'bitscore'], inplace=True)

    if df.empty:
        print("No valid E-value or bitscore data after conversion. No analysis can be performed.")
        return

    print(f"Number of loaded and valid hits: {len(df)}")

    # --- E-value Analysis ---
    print("\n--- E-value Analysis ---")

    # --- THIS IS THE MODIFIED PART ---
    # Find the minimum E-value in the entire dataset
    min_evalue = df['evalue'].min()
    print(f"Minimum E-value found in the null model: {min_evalue:.2e}")
    # ---------------------------------

    # E-value Histogram
    plt.figure(figsize=(12, 7)) # Increased figure size
    evalues_gt_zero = df['evalue'][df['evalue'] > 0]
    if not evalues_gt_zero.empty:
        min_val = max(1e-200, evalues_gt_zero.min())
        max_val = max(1.0, evalues_gt_zero.max())
        if min_val >= max_val:
            max_val = min_val * 10 if min_val > 0 else 1.0
            if min_val >= max_val:
                 min_val = max_val /100 if max_val > 0 else 1e-2

        bins = np.logspace(np.log10(min_val), np.log10(max_val), 50)
        
        plt.hist(evalues_gt_zero, bins=bins, color='skyblue', edgecolor='black', label='E-value Frequency')
        plt.xscale('log')
        plt.title('E-value Histogram of Null Model (E-values > 0, log scale on X-axis)', fontsize=14)

        # --- THIS IS THE MODIFIED PART ---
        # Add a vertical red dashed line at the minimum E-value
        plt.axvline(x=min_evalue, color='red', linestyle='--', linewidth=2, 
                    label=f'Minimum E-value ({min_evalue:.2e})')
        plt.legend() # Display the legend to show the line's label
        # ---------------------------------

    else:
        num_zeros = len(df[df['evalue'] == 0])
        plt.hist(df['evalue'], bins=10, color='skyblue', edgecolor='black')
        plt.title(f'E-value Histogram (all values; {num_zeros} are 0.0)', fontsize=14)

    plt.xlabel('E-value', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    plt.savefig('evalue_histogram.png')
    print("E-value histogram saved as 'evalue_histogram.png'")


    # --- Bitscore Analysis ---
    print("\n--- Bitscore Analysis ---")

    # Bitscore Histogram
    plt.figure(figsize=(12, 7))
    plt.hist(df['bitscore'], bins=50, color='lightgreen', edgecolor='black')
    plt.title('Bitscore Histogram of Null Model', fontsize=14)
    plt.xlabel('Bitscore', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    plt.savefig('bitscore_histogram.png')
    print("Bitscore histogram saved as 'bitscore_histogram.png'")

    # Calculate the maximum Bitscore
    max_bitscore = df['bitscore'].max()
    print(f"Maximum Bitscore found in the null model: {max_bitscore:.2f}")

if __name__ == '__main__':
    analyze_blast_output(file_path="uniprot_sprot_shuffle.blastp.out")
