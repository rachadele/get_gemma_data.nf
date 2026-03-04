# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A Nextflow pipeline for downloading and preprocessing single-cell RNA-seq datasets from the [GEMMA](https://gemma.msl.ubc.ca/) (staging) database. It downloads MEX expression matrices, cell type assignments, and sample metadata, then processes everything into AnnData (H5AD) format.

## Running the Pipeline

```bash
# Required env vars
export GEMMA_USERNAME="your_username"
export GEMMA_PASSWORD="your_password"

# Basic run from a study list file
nextflow run main.nf -profile conda --study_names study_names_human.txt

# Resume a failed/partial run
nextflow run main.nf -profile conda --study_names study_names_human.txt -resume

# Use a pre-downloaded studies directory instead of downloading
nextflow run main.nf -profile conda --study_paths /path/to/existing/studies

# Per-sample mode (one H5AD per sample instead of one per study)
nextflow run main.nf -profile conda --study_names study_names_human.txt --process_samples true

# Author-submitted cell types (default: false, uses curated assignments)
nextflow run main.nf -profile conda --study_names study_names_human.txt --author_submitted true
```

The `--outdir` is auto-generated from params: `{study_names}_author_{author_submitted}_process_samples_{process_samples}`.

## Key Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--study_names` | null | Text file, one GEMMA study ID per line |
| `--study_file` | null | Comma- or space-separated study IDs (inline) |
| `--study_paths` | null | Path to pre-downloaded MEX directory |
| `--process_samples` | false | true = one H5AD per sample |
| `--author_submitted` | false | Cell type assignment source |
| `--gene_mapping` | `meta/gemma_genes.tsv` | ENSEMBL_ID → OFFICIAL_SYMBOL |

Provide exactly one of `--study_names`, `--study_file`, or `--study_paths`.

## Architecture

```
main.nf                          # Main workflow entry point
nextflow.config                  # Params, SLURM executor config (queueSize=90)
modules/
  subworkflows/
    download_studies.nf          # Reads study list or discovers existing dirs
    process_queries.nf           # Unused subworkflow (superseded by main.nf logic)
  processes/
    download_studies.nf          # Runs gemma-cli-staging to fetch MEX data
    process_query_combined.nf    # PROCESS_QUERY_COMBINED: one H5AD per study
    process_query_samples.nf     # PROCESS_QUERY_SAMPLE: one H5AD per sample
bin/
  get_gemma_meta.py              # Fetches sample metadata via gemmapy API
  standardize_metadata.py        # Renames metadata columns to standard schema
  process_query.py               # Builds combined H5AD (called by PROCESS_QUERY_COMBINED)
  process_query_samples.py       # Builds per-sample H5AD (called by PROCESS_QUERY_SAMPLE)
  extract_downloaded_mex.sh      # Utility: copies MEX dirs from Nextflow work/ into a new dir
meta/
  gemma_genes.tsv                # Gene mapping: ENSEMBL_ID → OFFICIAL_SYMBOL
```

## Workflow Data Flow

1. **Download**: `DOWNLOAD_STUDIES_SUBWF` either reads study IDs from a file and runs `gemma-cli-staging`, or discovers existing directories from `--study-paths`. Emits `(study_name, study_dir)` tuples.
2. **Cell types**: `downloadCelltypes` fetches cell type assignments via GEMMA REST API (curl). The `author_submitted` param controls which protocol endpoint is used.
3. **Sample metadata**: `getGemmaMeta` calls `get_gemma_meta.py` (uses `gemmapy` Python client) to fetch and pivot sample characteristics.
4. **Metadata standardization**: `standardizeMetadata` renames columns using the `RENAME_MAP` in `standardize_metadata.py` (e.g., "biological sex" → "sex", "organism part" → "region").
5. **H5AD processing**: Depending on `process_samples`:
   - Combined: `PROCESS_QUERY_COMBINED` → `process_query.py` merges all samples per study
   - Sample: `PROCESS_QUERY_SAMPLE` → `process_query_samples.py` processes each sample independently; samples with <50 cells go to `small_samples/`

## Conda Environment

All Python processes use `/home/rschwartz/anaconda3/envs/scanpyenv`. Required packages: `scanpy`, `anndata`, `pandas`, `scipy`, `gemmapy`.

## Cluster Config

Runs on SLURM with `-C thrd64 --cpus-per-task=10`. Queue size is 90. Edit `nextflow.config` to change.

## Utility Scripts

`bin/extract_downloaded_mex.sh <work_dir> <new_dir>` — copies MEX directories out of Nextflow's `work/` directory into a flat structure. Useful for recovering data when re-running with `--study-paths`.

## Output Layout

```
{outdir}/
├── mex/                          # Raw MEX (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
├── cell_type_assignments/        # {study}.celltypes.tsv
├── metadata/{study}/{organism}/  # {study}_sample_meta.tsv
├── metadata_standardized/        # Column-renamed version of above
├── unique_cells/{study}/         # Cell type count summaries
├── h5ad/{study}/                 # {study}.h5ad or {study}_{sample}.h5ad
└── small_samples/{study}/        # Samples with <50 cells
```

## H5AD Schema

- **obs**: `sample_id`, `cell_type`, `cell_type_uri`, `region`, `sex`, `dev_stage`, `organism`, `assay`, `donor_id`, plus study-specific columns
- **var**: index = ENSEMBL_ID, `feature_name` = OFFICIAL_SYMBOL
- **X**: CSR sparse matrix (raw counts)
