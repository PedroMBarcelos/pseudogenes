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
- Input resolution from local files or download URLs.
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
./nextflow run main.nf -preview -profile local
```

3. Run locally:

```bash
./nextflow run main.nf -profile local -resume
```

4. Run with HPC profile (SLURM base settings):

```bash
./nextflow run main.nf -profile hpc -resume
```

5. Enable containers when desired:

```bash
./nextflow run main.nf -profile local,docker -resume
./nextflow run main.nf -profile hpc,singularity -resume
```

Important parameters (nextflow.config):
- params.uniprot_fasta
- params.genome_fasta
- params.blast_bin_dir
- params.ssearch_bin
- params.ssearch_chunk_size

Threshold source parameters:
- params.run_null_model (default: false)
- params.min_evalue_file (default: null)
- params.min_evalue_manual (default: null)

Threshold precedence:
1. If `min_evalue_file` is set, that file is used.
2. Else if `min_evalue_manual` is set, that value is used.
3. Else if `run_null_model` is true, blastp null-model is executed and min e-value is extracted.
4. Else the workflow exits with a configuration error.

Examples:

Run with default behavior (manual threshold required):

```bash
./nextflow run main.nf -profile local --min_evalue_manual 1e-20 -resume
```

Run null-model explicitly when desired:

```bash
./nextflow run main.nf -profile local --run_null_model true -resume
```

Skip null-model and use a manual threshold:

```bash
./nextflow run main.nf -profile local --run_null_model false --min_evalue_manual 1e-20 -resume
```

Skip null-model and reuse a previous threshold file:

```bash
./nextflow run main.nf -profile local --run_null_model false --min_evalue_file ./results_real/metrics/min_evalue.txt -resume
```

# How to run?
Use the Nextflow workflow instead of the legacy shell scripts.

```bash
./nextflow run main.nf -preview -profile local
./nextflow run main.nf -profile local -resume
```

The old shell entry points have been retired because their logic now lives in the DSL2 workflow.

