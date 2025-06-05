#!/bin/bash

#!/bin/bash

echo "Checking and downloading BLAST..."
if [ ! -f ncbi-blast-2.16.0+-x64-linux.tar.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
    tar -xzf ncbi-blast-2.16.0+-x64-linux.tar.gz
else
    echo "BLAST archive already exists. Skipping download."
fi

echo "Checking and downloading E. coli genome..."
if [ ! -f GCF_000005845.2_ASM584v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
else
    echo "E. coli genome already exists. Skipping download."
fi

echo "Checking and downloading UniProt..."
if [ ! -f uniprot_sprot.fasta.gz ]; then
    wget https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
else
    echo "UniProt file already exists. Skipping download."
fi

echo "Running Python script..."
python3 Python_translation.py uniprot_sprot.fasta.gz

echo "Running BLAST..."
# --- Configuration ---
BLAST_VERSION="2.16.0+" # Define the version you are using
BLAST_ARCHIVE_NAME="ncbi-blast-${BLAST_VERSION}-x64-linux.tar.gz"
BLAST_EXTRACTED_DIR="ncbi-blast-${BLAST_VERSION}"
MAKEBLASTDB_EXEC="${BLAST_EXTRACTED_DIR}/bin/makeblastdb"

# Input and output for makeblastdb (adjust these to your actual file names)
INPUT_FASTA="uniprot_sprot.fasta.gz.shuffle" 
DB_NAME="uniprot_sprot_shuffle"
DB_TYPE="prot"

# --- Check for BLAST+ and extract if necessary ---
echo "Checking for BLAST+ installation..."

if [ ! -f "${MAKEBLASTDB_EXEC}" ]; then
  echo "'${MAKEBLASTDB_EXEC}' not found."
  if [ -f "${BLAST_ARCHIVE_NAME}" ]; then
    echo "Found archive '${BLAST_ARCHIVE_NAME}'. Extracting..."
    tar -zxvf "${BLAST_ARCHIVE_NAME}"
    # Check if extraction was successful by looking for the executable again
    if [ -f "${MAKEBLASTDB_EXEC}" ]; then
      echo "BLAST+ extracted successfully."
    else
      echo "ERROR: Failed to find '${MAKEBLASTDB_EXEC}' after attempting extraction."
      echo "Please check the archive and extraction process."
      exit 1 # Exit the script if extraction failed
    fi
  else
    echo "ERROR: BLAST+ archive '${BLAST_ARCHIVE_NAME}' not found in the current directory."
    echo "Please download it or ensure it's in the correct location."
    exit 1 # Exit the script if archive is missing
  fi
else
  echo "BLAST+ executable '${MAKEBLASTDB_EXEC}' found. Skipping extraction."
fi

# --- Proceed with your pipeline steps ---
echo "Proceeding to create BLAST database..."

# Check if the input FASTA file exists
if [ ! -f "${INPUT_FASTA}" ]; then
    echo "ERROR: Input FASTA file '${INPUT_FASTA}' not found!"
    exit 1
fi

# Now, run makeblastdb using the variable for the executable
"${MAKEBLASTDB_EXEC}" -in "${INPUT_FASTA}" -dbtype "${DB_TYPE}" -out "${DB_NAME}"

if [ $? -eq 0 ]; then
  echo "BLAST database '${DB_NAME}' created successfully."
else
  echo "ERROR: makeblastdb command failed."
  exit 1
fi

echo "gunzip uniprot_sprot.fasta.gz..."

gunzip uniprot_sprot.fasta.gz

echo "Start bastp with nohup..."

ncbi-blast-2.16.0+/bin/blastp \
-num_threads 4 \
-word_size 3 \
-gapopen 11 \
-gapextend 1 \
-matrix BLOSUM62 \
-threshold 13 \
-comp_based_stats 2 \
-seg yes \
-soft_masking true \
-lcase_masking \
-evalue 10 \
-max_target_seqs 1000000 \
-outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovs qcovhsp" \
-query uniprot_sprot.fasta \
-db uniprot_sprot_shuffle \
-out uniprot_sprot_shuffle.blastp.out 


echo "Lowest E-Value in the Null Model..."
 awk '{print $13,$11}' uniprot_sprot_shuffle.blastp.out | LC_ALL=C sort -gr -k 2 | tail -n 1

echo "Calculating E-Value and Bitscore threshhold..."
python3 analise_ blast.py





