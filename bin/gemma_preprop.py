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

def load_mex(study_path):


  sample_ids = os.listdir(study_path)
  
  all_sample_ids = {}
  
  for sample_id in sample_ids:
    query_path = os.path.join(study_path, sample_id)
    new_sample_id = sample_id.split("_")[0]
    try:
        # Attempt to read the 10x mtx data
        adata = sc.read_10x_mtx(query_path)
        adata.obs_names_make_unique()
        all_sample_ids[new_sample_id] = adata
    except Exception as e:
        print(f"Error processing {sample_id}: {e}")
        
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
            adata.obs_names_make_unique()  # make sure the observation names are unique
            # Store the AnnData object in the dictionary
            all_sample_ids[new_sample_id] = adata
            print(f"Successfully created AnnData for {sample_id} from individual files.")
        
        except Exception as manual_e:
            print(f"Error processing {sample_id} manually: {manual_e}")
            all_sample_ids[new_sample_id] = None  # Or handle it differently, e.g., skip this sample
    
  return all_sample_ids


  
def combine_adata(all_sample_ids):
  combined_adata = sc.concat(all_sample_ids, label="sample_id", join="inner") 
  combined_adata.obs["cell_id"] = combined_adata.obs.index
  combined_adata.obs_names = combined_adata.obs["cell_id"].astype(str) + "_" + combined_adata.obs["sample_id"].astype(str)
  combined_adata.obs.index = combined_adata.obs_names
  #combined_adata.write_h5ad(f"{study_name}.h5ad")
  return(combined_adata)
  
def add_cell_meta(adata, cell_meta_path, sep="\t", outdir=None, study_name=None):
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
  #adata.write_h5ad(os.path.join(outdir, f"{study_name}.h5ad"))

def add_sample_meta(adata, meta_path, sep="\t", outdir=None, study_name=None):
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
  
def write_unique_cells(meta_path, sep='\t', outdir=None, study_name=None):
  meta = pd.read_csv(meta_path, sep=sep)
  unique_cells = pd.DataFrame(meta["cell_type"].value_counts()).reset_index()
  unique_cells.to_csv(f"{study_name}_unique_cells.tsv", sep='\t', index=False)

def write_as_samples(adata, study_name=None):
  adata.obs["sample_id"] = adata.obs["sample_id"].apply(str)
  for sample_id in adata.obs["sample_id"].unique():
    sample_adata = adata[adata.obs["sample_id"] == sample_id]
    #drop columns with only NaNs
    #fille NaN with ""
    sample_adata.obs = sample_adata.obs.fillna("")
    sample_adata.write_h5ad(f"{study_name}_{sample_id}.h5ad")
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study_dir", type=str, required=True)
    parser.add_argument("--study_name", type=str, required=True)
    parser.add_argument("--cell_meta_path", type=str, required=True)
    parser.add_argument("--sample_meta_path", type=str, required=True)
    parser.add_argument("--outdir", type=str, required=True)
    
    args = parser.parse_args()
    study_path = args.study_dir
    study_name = args.study_name
    cell_meta_path = args.cell_meta_path
    sample_meta_path = args.sample_meta_path
    outdir = args.outdir
    
    all_sample_ids = load_mex(study_path)
    adata = combine_adata(all_sample_ids)
    adata = add_cell_meta(adata, cell_meta_path, study_name=study_name, outdir=outdir)
    adata = add_sample_meta(adata, sample_meta_path, study_name=study_name, outdir=outdir)
    write_as_samples(adata, study_name=study_name)

    write_unique_cells(cell_meta_path, study_name=study_name)
    
if __name__ == "__main__":
    main()
    