#!/bin/bash

echo "Downloading BLAST..."
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz

echo "Downloading E. coli genome..."
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/005/845/GCF_037005845.1_ASM3700584v1/GCF_037005845.1_ASM3700584v1_genomic.fna.gz

echo "Downloading UniProt..."
wget https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

echo "Running Python script..."
python3 Python_translation.py uniprot_sprot.fasta.gz
