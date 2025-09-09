#!/bin/bash

# Função para verificar se o arquivo existe e perguntar se deseja sobrescrever
check_overwrite() {
    FILE_TO_CHECK=$1
    if [ -f "$FILE_TO_CHECK" ]; then
        # Aguarda 3 minutos (180s) por uma resposta. Se não houver, continua.
        read -t 180 -p "O arquivo '$FILE_TO_CHECK' já existe. Deseja sobrescrevê-lo? (y/n) [Padrão: y após 3 min] " -n 1 -r
        read_exit_status=$?
        echo

        # Se o read estourou o tempo (exit status > 128) ou o usuário apenas apertou Enter, prossiga.
        if [ $read_exit_status -ne 0 ]; then
            echo "Nenhuma resposta recebida. Prosseguindo com a sobrescrita."
            return 0 # Retorna 0 para indicar que pode continuar
        fi

        # Se o usuário respondeu, verifique a resposta.
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Operação cancelada pelo usuário. Pulando esta etapa."
            return 1 # Retorna 1 para indicar que não deve continuar
        fi
    fi
    return 0 # Retorna 0 para indicar que pode continuar
}


# Pastas organizadas
BLAST_DIR="BLAST"
GENOME_DIR="genoma"
RESULTS_DIR="resultados"

# Set NUM_JOBS to be the total number of cores minus one.
total_cores=$(nproc)
if (( total_cores > 1 )); then
    NUM_JOBS=$((total_cores - 1))
else
    NUM_JOBS=1
fi
echo "O numeros de threads disponíveis é "$NUM_JOBS""

mkdir -p "$BLAST_DIR" "$GENOME_DIR" "$RESULTS_DIR"

echo "Checking and downloading BLAST..."
if [ ! -f "$BLAST_DIR/ncbi-blast-2.17.0+-x64-linux.tar.gz" ]; then
    wget -O "$BLAST_DIR/ncbi-blast-2.17.0+-x64-linux.tar.gz" https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
    tar -xzf "$BLAST_DIR/ncbi-blast-2.17.0+-x64-linux.tar.gz" -C "$BLAST_DIR"
else
    echo "BLAST archive already exists. Skipping download."
fi

echo "Checking and downloading E. coli genome..."
if [ ! -f "$GENOME_DIR/GCF_000005845.2_ASM584v2_genomic.fna.gz" ]; then
    wget -O "$GENOME_DIR/GCF_000005845.2_ASM584v2_genomic.fna.gz" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
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

# Configuração dos caminhos dos executáveis
BLAST_VERSION="2.17.0+"
BLAST_EXTRACTED_DIR="$BLAST_DIR/ncbi-blast-${BLAST_VERSION}"
MAKEBLASTDB_EXEC="${BLAST_EXTRACTED_DIR}/bin/makeblastdb"
BLASTP_EXEC="${BLAST_EXTRACTED_DIR}/bin/blastp"
TBLASTN_EXEC="${BLAST_EXTRACTED_DIR}/bin/tblastn"

# Input e output para makeblastdb
INPUT_FASTA="uniprot_sprot.fasta.gz.shuffle"
DB_NAME="${RESULTS_DIR}/uniprot_sprot_shuffle"
DB_TYPE="prot"

echo "Checking for BLAST+ installation..."
if [ ! -f "${MAKEBLASTDB_EXEC}" ]; then
    echo "'${MAKEBLASTDB_EXEC}' not found."
    if [ -f "$BLAST_DIR/ncbi-blast-${BLAST_VERSION}-x64-linux.tar.gz" ]; then
        echo "Found archive. Extracting..."
        tar -zxvf "$BLAST_DIR/ncbi-blast-${BLAST_VERSION}-x64-linux.tar.gz" -C "$BLAST_DIR"
        if [ ! -f "${MAKEBLASTDB_EXEC}" ]; then
            echo "ERROR: Failed to find '${MAKEBLASTDB_EXEC}' after extraction."
            exit 1
        fi
    else
        echo "ERROR: BLAST+ archive not found."
        exit 1
    fi
else
    echo "BLAST+ executable found. Skipping extraction."
fi

echo "Proceeding to create BLAST database..."
if [ ! -f "${INPUT_FASTA}" ]; then
    echo "ERROR: Input FASTA file '${INPUT_FASTA}' not found!"
    exit 1
fi

# Medir o tempo de criação do banco de dados UniProt
start_time=$(date +%s)

if ! check_overwrite "${DB_NAME}.psq"; then
    echo "Criação do banco de dados UniProt pulada."
else
    "${MAKEBLASTDB_EXEC}" -in "${INPUT_FASTA}" -dbtype "${DB_TYPE}" -out "${DB_NAME}"
    if [ $? -eq 0 ]; then
      echo "BLAST database '${DB_NAME}' created successfully."
    else
      echo "ERROR: makeblastdb command failed."
      exit 1
    fi
fi

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Criação do banco de dados UniProt concluída em $duration segundos."

echo "gunzip uniprot_sprot.fasta.gz..."
gunzip -f uniprot_sprot.fasta.gz

echo "Start blastp..."
# Medir o tempo do blastp
start_time=$(date +%s)

if ! check_overwrite "${RESULTS_DIR}/uniprot_sprot_shuffle.blastp.out"; then
    echo "Execução do BLASTP pulada."
else
    "${BLASTP_EXEC}" \
    -num_threads $NUM_JOBS \
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
    -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovs qcovhsp" \
    -query uniprot_sprot.fasta \
    -db "${DB_NAME}" \
    -out "${RESULTS_DIR}/uniprot_sprot_shuffle.blastp.out"
fi

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "BLASTP concluído em $duration segundos."


echo "Genome-Wide Homology-Based search..."

if [ ! -f "$GENOME_DIR/GCF_000005845.2_ASM584v2_genomic.fna" ]; then
    echo "Unzipping E. coli genome..."
    gunzip -k "$GENOME_DIR/GCF_000005845.2_ASM584v2_genomic.fna.gz"
fi

# Medir o tempo de criação do banco de dados do genoma
start_time=$(date +%s)

if ! check_overwrite "${RESULTS_DIR}/ecolik12db.nsq"; then
    echo "Criação do banco de dados do genoma pulada."
else
    "${MAKEBLASTDB_EXEC}" -in "$GENOME_DIR/GCF_000005845.2_ASM584v2_genomic.fna" -dbtype nucl -out "${RESULTS_DIR}/ecolik12db"
fi

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Criação do banco de dados do genoma concluída em $duration segundos."

# Medir o tempo do tblastn
start_time=$(date +%s)

if ! check_overwrite "${RESULTS_DIR}/genome.tblastn.out"; then
    echo "Execução do TBLASTN pulada."
else
    "${TBLASTN_EXEC}" \
    -num_threads $NUM_JOBS \
    -word_size 3 \
    -gapopen 11 \
    -gapextend 1 \
    -matrix BLOSUM62 \
    -threshold 13 \
    -comp_based_stats 2 \
    -seg yes \
    -soft_masking true \
    -lcase_masking \
    -evalue 1 \
    -max_target_seqs 1000000 \
    -dbsize 10000 \
    -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sframe sstrand qcovs qcovhsp" \
    -query uniprot_sprot.fasta \
    -db "${RESULTS_DIR}/ecolik12db" \
    -out "${RESULTS_DIR}/genome.tblastn.out"
fi

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "TBLASTN concluído em $duration segundos."