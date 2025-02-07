import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from pathlib import Path
import argparse
import os
import json
import sys
import numpy as np
import gemmapy

def argument_parser():
    parser = argparse.ArgumentParser(description="Preprocess data from GEMMA")
    parser.add_argument("--study_name", type=str, help="Name of the study")
    return parser.parse_args()

def main():
    args = argument_parser()
   # username = os.getenv("GEMMA_USERNAME")
    #password = os.getenv("GEMMA_PASSWORD")
    client = gemmapy.GemmaPy(auth=["raschwar","7nddtt"], path='dev')
    study_name = args.study_name
    samples_raw = client.raw.get_dataset_samples(study_name)
  
    samples = client.get_dataset_samples(study_name)
    sample_ids = [x.id for x in samples_raw.data]
    organisms = [x.array_design.taxon.scientific_name.lower().replace(" ", "_") for x in samples_raw.data]
    sample_meta = [df for df in samples["sample_characteristics"]]
    # add sample id to each df in sample meta
    sample_meta_updated = [df.assign(sample_id=sample_ids[i]) for i, df in enumerate(sample_meta)]
    # add organism to each df in sample meta
    sample_meta_updated = [df.assign(organism=organisms[i]) for i, df in enumerate(sample_meta_updated)]
    # combine all dfs
    sample_meta_combined = pd.concat(sample_meta_updated)
    sample_meta_combined.drop_duplicates(subset=["sample_id", "category"], inplace=True)

    sample_meta_df = sample_meta_combined.pivot(index=["sample_id"], columns="category",values="value").reset_index()
    sample_meta_df.to_csv(f"{study_name}_sample_meta.tsv", index=False, sep="\t")
    

    
if __name__ == "__main__":
    main()
