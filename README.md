# pseudogenes

This repository contains a workflow intended to identify and analyze pseudogenes from a givem genome sequence.

The pipeline.md file is the original code from Marcos, it has been updated to get functional links.
The code in perl has been translated to python in the Python_translation.py.

## Nextflow implementation (DSL2)

A first Nextflow DSL2 implementation is now available.

Main files:
- main.nf
- nextflow.config
- workflows/pseudogenes.nf
- modules/local/*.nf

What is already implemented:
- Input resolution for local files (genome input is required from user path).
- UniProt preparation and shuffled sequence generation.
- BLAST database creation (protein and nucleotide).
- Optional blastp NULL model and dynamic minimum e-value extraction.
- tblastn genome search and threshold-based filtering.
- SSEARCH fragment chunking, scatter processing, and gather.
- BED conversion, sorting, merge, and final report generation.

Quick start:
1. Ensure Nextflow is available (or use ./nextflow if present in this folder).
2. Run preview mode (no heavy execution):

```bash
./nextflow run main.nf -preview -profile local --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna
```

3. Run locally:

```bash
./nextflow run main.nf -profile local --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna -resume
```

4. Run with HPC profile (SLURM base settings):

```bash
./nextflow run main.nf -profile hpc --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna -resume
```

5. Enable containers when desired:

```bash
./nextflow run main.nf -profile local,docker --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna -resume
./nextflow run main.nf -profile hpc,singularity --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna -resume
```

Important parameters (nextflow.config):
- params.uniprot_fasta
- params.genome_fasta
- params.blast_bin_dir
- params.ssearch_bin
- params.ssearch_chunk_size

Threshold source parameters:
- params.min_evalue (default: 1e-7)
- params.run_null_model (default: false)

Genome input parameter:
- params.genome_fasta (required; local FASTA or gzipped FASTA)

Threshold behavior:
1. By default, `min_evalue` is used (default `1e-7`).
2. If `run_null_model` is true, blastp null-model is executed and the generated minimum e-value is used.

Genome behavior:
1. `genome_fasta` must be provided by the user.
2. If the input ends with `.gz`, it is decompressed automatically.
3. If not compressed, it is used directly.

Examples:

Run with default threshold (`1e-7`):

```bash
./nextflow run main.nf -profile local --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna -resume
```

Run with custom threshold:

```bash
./nextflow run main.nf -profile local --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna --min_evalue 1e-20 -resume
```

Run null-model explicitly (uses generated threshold):

```bash
./nextflow run main.nf -profile local --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna --run_null_model true -resume
```

Run with gzipped genome input:

```bash
./nextflow run main.nf -profile local --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna.gz --min_evalue 1e-20 -resume
```

# How to run?
Use the Nextflow workflow instead of the legacy shell scripts.

```bash
./nextflow run main.nf -preview -profile local --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna
./nextflow run main.nf -profile local --genome_fasta ./genoma/GCF_000005845.2_ASM584v2_genomic.fna -resume
```

The old shell entry points have been retired because their logic now lives in the DSL2 workflow.

