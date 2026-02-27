# GEMMA Single-Cell Data Download Pipeline

A Nextflow pipeline for downloading and preprocessing single-cell RNA-seq datasets from the [GEMMA](https://gemma.msl.ubc.ca/) database. The pipeline downloads expression matrices, cell type assignments, and sample metadata, then processes them into standardized AnnData (H5AD) format.

## Features

- Downloads single-cell expression data in MEX format from GEMMA
- Retrieves cell type assignments and sample metadata via GEMMA API
- Standardizes gene names using ENSEMBL ID mapping
- Supports two processing modes:
  - **Combined mode**: All samples merged into a single H5AD per study
  - **Sample mode**: Each sample processed into a separate H5AD file
- SLURM cluster support with conda environment management

## Requirements

### Software
- [Nextflow](https://www.nextflow.io/) (>=21.04)
- [Conda](https://docs.conda.io/)
- `gemma-cli-staging` (GEMMA command-line tool)
- `curl`

### Python Dependencies
The pipeline requires a conda environment with:
- scanpy
- pandas
- anndata
- scipy
- gemmapy

### Credentials
GEMMA API credentials must be set as environment variables:
```bash
export GEMMA_USERNAME="your_username"
export GEMMA_PASSWORD="your_password"
```

## Input

### Study Names File
A text file with one GEMMA study ID per line:
```
GSE237718
GSE180670
Velmeshev-2019.1
```

Example files are provided:
- `study_names_human.txt` - Human studies
- `study_names_mouse.txt` - Mouse studies

### Gene Mapping File
Located at `meta/gemma_genes.tsv`, containing ENSEMBL_ID to OFFICIAL_SYMBOL mappings.

## Usage

### Basic Usage

```bash
# Using a study names file
nextflow run main.nf -profile conda \
  --study-names study_names_human.txt

# Using an existing studies directory
nextflow run main.nf -profile conda \
  --study-paths /path/to/existing/studies
```

### Processing Modes

```bash
# Combined mode (default) - one H5AD per study
nextflow run main.nf -profile conda \
  --study-names study_names_human.txt \
  --process_samples false

# Sample mode - one H5AD per sample
nextflow run main.nf -profile conda \
  --study-names study_names_human.txt \
  --process_samples true
```

### Additional Options

```bash
nextflow run main.nf -profile conda \
  --study-names study_names_human.txt \
  --author_submitted true \
  --outdir custom_output_dir \
  -resume
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--study-names` | Path to file with study IDs (one per line) | `null` |
| `--study-file` | Comma- or space-separated study IDs (inline) | `null` |
| `--study-paths` | Path to pre-downloaded studies directory | `null` |
| `--process_samples` | Process each sample separately (`true`) or combined (`false`) | `false` |
| `--author_submitted` | Use author-submitted cell types | `false` |
| `--gene_mapping` | Path to gene mapping TSV | `meta/gemma_genes.tsv` |
| `--outdir` | Output directory | Auto-generated based on params |

**Note:** You must provide exactly one of `--study-names`, `--study-file`, or `--study-paths`.

## Output Structure

```
{outdir}/
├── mex/                          # Raw MEX format data
│   └── {study}/
│       └── {sample}/
│           ├── matrix.mtx.gz
│           ├── features.tsv.gz
│           └── barcodes.tsv.gz
├── cell_type_assignments/        # Cell type annotations
│   └── {study}.celltypes.tsv
├── metadata/                     # Sample metadata
│   └── {study}/
│       └── {organism}/
│           └── {study}_sample_meta.tsv
├── unique_cells/                 # Cell type count summaries
│   └── {study}/
│       └── {study}_unique_cells.tsv
├── h5ad/                         # Processed AnnData files
│   └── {study}/
│       └── {study}.h5ad          # or {study}_{sample}.h5ad
└── small_samples/                # Samples with <50 cells
    └── {study}/
```

## Output File Formats

### H5AD Files
AnnData format containing:
- **X**: Sparse expression matrix (CSR format)
- **obs**: Cell metadata (sample_id, cell_type, region, sex, dev_stage, organism, assay)
- **var**: Gene annotations (feature_name, ENSEMBL_ID)

### Cell Type Assignments TSV
| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| cell_id | Cell barcode |
| cell_type | Assigned cell type |
| cell_type_uri | Ontology URI |

### Sample Metadata TSV
Contains sample characteristics and factor values including:
- organism_part / region
- sex
- developmental_stage
- disease
- assay

## Pipeline Workflow

```
1. Input Handling
   └── Read study names OR discover from existing directory

2. Data Download
   ├── Download MEX expression matrices (gemma-cli-staging)
   ├── Download cell type assignments (GEMMA API)
   └── Extract sample metadata (gemmapy)

3. Processing
   ├── Combined mode: Merge all samples per study
   └── Sample mode: Process each sample separately

4. Output Generation
   ├── Standardize gene names
   ├── Integrate cell/sample metadata
   └── Export to H5AD format
```

## Cluster Configuration

The pipeline is configured for SLURM execution:
- Queue size: 90 concurrent jobs
- CPUs per task: 10
- Cluster options: `-C thrd64 --cpus-per-task=10`

To modify, edit `nextflow.config`:
```groovy
process {
    executor = 'slurm'
    queueSize = 90
    clusterOptions = '-C thrd64 --cpus-per-task=10'
}
```

## Resuming Failed Runs

Use the `-resume` flag to continue from the last successful checkpoint:
```bash
nextflow run main.nf -profile conda --study-names study_names.txt -resume
```

## Troubleshooting

### Authentication Errors
Ensure `GEMMA_USERNAME` and `GEMMA_PASSWORD` environment variables are set correctly.

### Small Samples
Samples with fewer than 50 cells are automatically moved to the `small_samples/` directory.

### Missing Gene Mappings
Genes without ENSEMBL ID mappings will retain their original identifiers.

## License

This project is developed for research purposes at UBC.

## Contact

For questions about GEMMA data, visit: https://gemma.msl.ubc.ca/
