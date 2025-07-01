#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from pathlib import Path
import argparse
import os
import json
import scipy
import gzip
    
def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("--query_path", type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/work/1d/93b79a110ad273088caf132ec2a4ac/1293348_H21.33.038")
  parser.add_argument("--query_name", type=str, default="1293348_H21.33.038")
  parser.add_argument("--cell_meta_path", type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/work/1d/93b79a110ad273088caf132ec2a4ac/SEA-AD-DLPFC-2024.celltypes.tsv")
  parser.add_argument("--sample_meta_path", type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/work/1d/93b79a110ad273088caf132ec2a4ac/SEA-AD-DLPFC-2024_sample_meta.tsv")
 # parser.add_argument("--write_samples", action="store_true", help="Write samples as individual files")
  parser.add_argument('--gene_mapping', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv", help='Path to the gene mapping file')  
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
     #   all_sample_ids[new_sample_id] = adata
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
      # write a fake h5ad to trick nextflow
      adata.write_h5ad(f"{sample_id}_empty.h5ad")
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
  
def write_as_samples(adata, query_name, organism):
  small_samples = []
  
  for sample_id in adata.obs["sample_id"].unique():
    sample_adata = adata[adata.obs["sample_id"] == sample_id]
    print(sample_adata.shape)
    # check size
    if not check_size(sample_adata):
      # nee to combine with other samples
      small_samples.extend(sample_id)
      continue
    #drop columns with only NaNs
    #fille NaN with ""
    sample_adata.obs = sample_adata.obs.fillna("")
    sample_adata.write_h5ad(f"{query_name}_{sample_id}.h5ad")
    
  if len(small_samples) > 0:
    print(f"Samples {small_samples} are too small to be written as individual files. They will be combined with other samples.")
    # combine small samples
    small_samples_adata = adata[adata.obs["sample_id"].isin(small_samples)]
    small_samples_adata.obs = small_samples_adata.obs.fillna("")
    if check_size(small_samples_adata):
      small_samples_adata.write_h5ad(f"{query_name}_small_samples.h5ad")
    else:
      print(f"Combined small samples are still too small to be written as a single file. Skipping writing.")
    

def map_genes(query, gene_mapping):
    # Check if the "feature_name" column exists in query.var
  if "feature_name" not in query.var.columns:
      # Merge gene_mapping with query.var based on the index
      query.var = query.var.merge(gene_mapping["OFFICIAL_SYMBOL"], left_index=True, right_index=True, how="left")
      # Rename the merged column to "feature_name"
      query.var.rename(columns={"OFFICIAL_SYMBOL": "feature_name"}, inplace=True)
  query.var_names = query.var.index
  return query

def main():

    args = parse_args()
    #query_path = args.study_dir
    query_path = args.query_path
    query_name = args.query_name
    cell_meta_path = args.cell_meta_path
    sample_meta_path = args.sample_meta_path
   # write_samples = args.write_samples
    gene_mapping = args.gene_mapping
    gene_mapping = pd.read_csv(gene_mapping, sep="\t", header=0)
    # Drop rows with missing values in the relevant columns
    gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
    # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
    gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
    gene_mapping.set_index("ENSEMBL_ID", inplace=True)

    
  #  all_sample_ids = load_mex(query_path)
    sample_id = query_name.split("_")[0]  # Extract the sample ID from the query name
    adata = load_mex(query_path, sample_id)
    #adata = combine_adata(all_sample_ids)
    adata = add_cell_meta(adata, cell_meta_path)
    adata = add_sample_meta(adata, sample_meta_path)
    adata = map_genes(adata, gene_mapping)
    if check_size(adata) is False:
      os.makedirs("small_samples", exist_ok=True)
      adata.write_h5ad(os.path.join("small_samples",f"{query_name}.h5ad"))
    else:
      adata.write_h5ad(f"{query_name}.h5ad")
        
if __name__ == "__main__":
    main()
    