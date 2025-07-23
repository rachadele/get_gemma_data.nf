#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
from scipy.sparse import csr_matrix
import warnings
import cellxgene_census
import cellxgene_census.experimental
import scvi
from pathlib import Path
import argparse
import os
import json
import scipy
import gzip
import random
import numpy as np

def parse_arguments():
  parser = argparse.ArgumentParser()
  parser.add_argument("--study_path", type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/work/1d/93b79a110ad273088caf132ec2a4ac/1293348_H21.33.038")
  parser.add_argument("--cell_meta_path", type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/work/1d/93b79a110ad273088caf132ec2a4ac/SEA-AD-DLPFC-2024.celltypes.tsv")
  parser.add_argument("--sample_meta_path", type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/work/1d/93b79a110ad273088caf132ec2a4ac/SEA-AD-DLPFC-2024_sample_meta.tsv")
  parser.add_argument('--gene_mapping', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv", help='Path to the gene mapping file')  
  parser.add_argument("--study_name", type=str, default="SEA-AD-DLPFC-2024", help="Name of the study for output files")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args


def load_mex(query_path, sample_id):
  
  try:
    # Attempt to read the 10x mtx data
    adata = sc.read_10x_mtx(query_path)
    adata.obs["sample_id"] = sample_id  # Add sample_id to obs
    adata.obs.index = adata.obs_names + "_" + adata.obs["sample_id"]
    adata.obs_names_make_unique()
  except Exception as e:
    print(f"Error processing {sample_id} automatically: {e}. Trying manual read.")
    
    # If an error occurs, try reading the files manually
    try:
        # Read the matrix, genes, and barcodes files manually
        matrix_path = os.path.join(query_path, "matrix.mtx.gz")
        genes_path = os.path.join(query_path, "features.tsv.gz")
        barcodes_path = os.path.join(query_path, "barcodes.tsv.gz")
        
        # Load the matrix in CSR format
        with gzip.open(matrix_path, 'rb') as f:
            matrix = scipy.io.mmread(f).tocsr()
        
        # Read the gene and barcode files
        with gzip.open(genes_path, 'rt') as f:
            genes = [line.strip().split("\t") for line in f]
        with gzip.open(barcodes_path, 'rt') as f:
            barcodes = [line.strip() for line in f]
        
        # Create AnnData object
        adata = sc.AnnData(X=matrix.T)  # Transpose to match expected shape (cells x genes)
        adata.var_names = [gene[1] for gene in genes]
        adata.var_names_make_unique()# gene ids as the variable names
        adata.obs_names = barcodes  # cell barcodes as the observation names
        adata.obs["sample_id"] = sample_id  # Add sample_id to obs
        adata.obs.index = adata.obs_names + "_" + adata.obs["sample_id"]
        adata.obs_names_make_unique()  # make sure the observation names are unique
        # Store the AnnData object in the dictionary
     #   all_samples[new_sample_id] = adata
        print(f"Successfully created AnnData for {sample_id} from individual files.")
    
    except Exception as manual_e:
      adata = sc.AnnData(X=csr_matrix((0, 0)))
      adata.obs["sample_id"] = sample_id  # Add sample_id to obs
      adata.obs.index = adata.obs_names + "_" + adata.obs["sample_id"]
      adata.write_h5ad(f"{sample_id}_empty.h5ad")
      #print(f"Error processing {sample_id} manually: {manual_e}")
      raise Warning(f"Failed to process {sample_id} using both automatic and manual methods: {manual_e}")

    if has_expression_data(adata) is False:
      adata = sc.AnnData(X=csr_matrix((0, 0)))
      # create a fake h5ad to trick
      raise Warning(f"Sample {sample_id} has no expression data. Skipping.")
  return adata

def has_expression_data(adata):
    return hasattr(adata, "X") and getattr(adata.X, "nnz", None) not in (0, None)

def check_size(adata):
  if adata.shape[0] < 50:
    return False
  else:
    return True
  
  
def add_cell_meta(adata, cell_meta_path, sep="\t"):
  meta = pd.read_csv(cell_meta_path, sep=sep)
  meta["combined_id"] = meta["cell_id"].astype(str) + "_" + meta["sample_id"].astype(str)
  meta = meta.set_index("combined_id")
  # Filter `meta` to exclude columns that overlap with `adata.obs`
  meta_cleaned = meta.loc[:, ~meta.columns.isin(adata.obs.columns)]
  # Perform the join operation
  adata.obs = adata.obs.join(meta_cleaned, how="left")
  # drop NaNs in cell_type column
  adata = adata[~adata.obs["cell_type"].isna()]
  return(adata)

def add_sample_meta(adata, meta_path, sep="\t"):
  meta = pd.read_csv(meta_path, sep=sep)
  meta["sample_id"] = meta["sample_id"].apply(str)
  # Perform the join operation
  adata.obs = adata.obs.merge(meta, left_on="sample_id", right_on="sample_id", how="left", suffixes=("", "_y"))
  adata.obs.rename(columns={"organism part": "region", "developmental_stage": "dev_stage", "biological sex": "sex"}, inplace=True)
  try:
    adata.obs["region"]=adata.obs["region"].replace("adult ","")
  except:
    pass
  return(adata)


def map_genes(query, gene_mapping):
  # Drop rows with missing values in the relevant columns
  gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
  # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
  gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
  gene_mapping.set_index("ENSEMBL_ID", inplace=True)

    # Check if the "feature_name" column exists in query.var
  if "feature_name" not in query.var.columns:
      # Merge gene_mapping with query.var based on the index
      query.var = query.var.merge(gene_mapping["OFFICIAL_SYMBOL"], left_index=True, right_index=True, how="left")
      # Rename the merged column to "feature_name"
      query.var.rename(columns={"OFFICIAL_SYMBOL": "feature_name"}, inplace=True)
  query.var_names = query.var.index
  return query



def main():

  # Parse command line arguments
  args = parse_arguments()
  #SEED = args.seed
  #random.seed(SEED)         # For `random`
  #np.random.seed(SEED)      # For `numpy`
  #scvi.settings.seed = SEED # For `scvi`
  study_path = args.study_path
  study_name = args.study_name
  cell_meta_path = args.cell_meta_path
  sample_meta_path = args.sample_meta_path
  gene_mapping_path = args.gene_mapping
  gene_mapping = pd.read_csv(gene_mapping_path, sep="\t", header=0)
  # Set organism and census_version from arguments
  sample_ids = os.listdir(study_path)
  
  all_samples = {}
  
  for sample_id in sample_ids:
    query_path = os.path.join(study_path, sample_id)
    new_sample_id = sample_id.split("_")[0]
    all_samples[new_sample_id] = load_mex(query_path, new_sample_id)
  
  combined_adata = sc.concat(all_samples, label="sample_id", join="inner") 
  combined_adata.obs["cell_id"] = combined_adata.obs.index
  combined_adata.obs_names = combined_adata.obs["cell_id"].astype(str) + "_" + combined_adata.obs["sample_id"].astype(str)
  # save unprocessed adata

  combined_adata = map_genes(combined_adata, gene_mapping)
  
  # Add cell metadata
  combined_adata = add_cell_meta(combined_adata, cell_meta_path)
  
  # Add sample metadata
  combined_adata = add_sample_meta(combined_adata, sample_meta_path)
  combined_adata.obs = combined_adata.obs.astype(str) 
  
  combined_adata.write_h5ad(f"{study_name}.h5ad")

  
  
if __name__ == "__main__":
    main()