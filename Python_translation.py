# DATA SOURCE

## blast-3.15.0+
    #wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-3.16.0+-x64-linux.tar.gz

## Escherichia coli str. K-13 substr. MG1655 reference genome
    #wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037777777777/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

## UniProtKB/Swiss-Prot 2024_02 of UniProtKB, published on Wed Apr 23 2025 (573,230 sequence entries)
    #wget https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz



import sys  # For reading command-line arguments (like input file name)
import random  # For shuffling the artificial sequence
from collections import defaultdict  # Provides a dictionary with default values (avoids key errors)
import gzip
# Function to read sequences from a multi-FASTA file
def parse_fasta(file_path):
    sequences = []  # List to store sequences from the file
    with gzip.open(file_path, "rt") as f:
        current_seq = ''  # Initialize a string to hold the current sequence
        for line in f:  # Read file line by line
            line = line.strip()  # Remove newline characters and extra spaces
            if line.startswith('>'):  # Header line (sequence ID or description)
                if current_seq:  # If there's a sequence already built
                    sequences.append(current_seq)  # Add it to the list
                    current_seq = ''  # Reset for next sequence
            elif line:  # If it's a sequence line (not empty)
                current_seq += line  # Append sequence characters
        if current_seq:  # Add the last sequence after file ends
            sequences.append(current_seq)
    return sequences  # Return list of sequences

# Function to count each amino acid and the total sequence length
def count_amino_acids(sequences):
    aa_counts = defaultdict(int)  # Dictionary to count each amino acid (default 0)
    total_length = 0  # Sum of all sequence lengths
    for seq in sequences:  # Go through each sequence
        for aa in seq:  # Go through each amino acid
            aa_counts[aa] += 1  # Count the amino acid
        total_length += len(seq)  # Add sequence length to total
    return aa_counts, total_length  # Return dictionary and total length

#number_of_sequences = len(sequences) # Get number of sequences 
#avg_seq_length = total_length / number_of_sequences # Calculate the average seq lenght by dividing total length to number of sequences


# Function to calculate average frequency of each amino acid
# and how many times each should appear in the artificial sequence
def calculate_average_aa_composition(aa_counts, total_length, avg_seq_length):
    avg_freq = {}  # Dictionary to store frequency of each amino acid
    aa_numbers = {}  # Dictionary to store how many times each should appear
    for aa, count in aa_counts.items():  # For each amino acid
        freq = count / total_length  # Frequency = count / total number of amino acids
        avg_freq[aa] = freq  # Store frequency
        aa_numbers[aa] = round(freq * avg_seq_length)  # Estimate how many should appear in average-length sequence
    return avg_freq, aa_numbers  # Return both dictionaries

# Function to build a list of amino acids forming one artificial sequence
def create_artificial_sequence(aa_numbers):
    artificial_seq = []  # List to hold the artificial sequence
    for aa, count in aa_numbers.items():  # For each amino acid
        artificial_seq.extend([aa] * count)  # Add 'count' copies of the amino acid
    return artificial_seq  # Return the list (sequence)

# Main function that coordinates the entire process
def main(infile):
    sequences = parse_fasta(infile)  # Read sequences from the input FASTA file
    number_of_seqs = len(sequences)  # Count how many sequences were read

    if number_of_seqs == 0:  # If no sequences found, stop the program
        print("No sequences found in the input file.")
        return

    aa_counts, total_seq_size = count_amino_acids(sequences)  # Count AAs and total size of all sequences

    average_seq_size = total_seq_size / number_of_seqs  # Calculate average sequence length

    # Get the frequency of each AA and how many should appear in an average-length artificial sequence
    avg_freq, aa_numbers = calculate_average_aa_composition(aa_counts, total_seq_size, average_seq_size)

    # Print some summary statistics
    print(f"\nNUMBER OF SEQUENCES: {number_of_seqs}")
    print(f"\nNUMBER OF AMINO ACIDS: {total_seq_size}\n")

    print("AVERAGE AMINO ACID FREQUENCY:")
    for aa, freq in avg_freq.items():  # Print each amino acid and its average frequency
        print(f"{aa} -> {freq}")
        
    print(f"\nAVERAGE SEQUENCE SIZE: {average_seq_size:.2f}\n")

    print("NUMBER OF AMINO ACIDS IN ARTIFICIAL SEQUENCE:")
    for aa, count in aa_numbers.items():  # Print each amino acid and its estimated number in an artificial sequence
        print(f"{aa} -> {count}")

    artificial_seq = create_artificial_sequence(aa_numbers)  # Create one artificial average sequence

    # Print the unshuffled artificial sequence
    print("\nEXAMPLE ARTIFICIAL SEQUENCE:\n")
    print("".join(artificial_seq))  # Convert list to string

    # Open output file for writing the shuffled sequences
    with open(f"{infile}.shuffle", 'w') as out:
        for i in range(1, 10001):  # Create 10,000 shuffled versions of the sequence
            random_seq = random.sample(artificial_seq, len(artificial_seq))  # Shuffle the sequence
            out.write(f">shuffled_{i}\n{''.join(random_seq)}\n")  # Write to output file in FASTA format



if __name__ == "__main__":
    if len(sys.argv) != 2:  # Check if user gave exactly one file name
        print("Usage: python shuffle.py <Multi-FASTA file>")  # Print usage instruction
        sys.exit(1)  # Exit with error
    main(sys.argv[1])  # Run the main function with the file path



