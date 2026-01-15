#!/bin/python

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
    parser.add_argument("--study_name", type=str, help="Name of the study", default="GSE181021.2")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
    
    
    
def extract_meta(df_list, sample_ids, sample_names):
    values = [df.assign(sample_id=sample_ids[i]) for i, df in enumerate(df_list)]
    #factor_values = [df.assign(organism=organisms[i]) for i, df in enumerate(factor_values)]
    values = [df.assign(sample_name=sample_names[i]) for i, df in enumerate(values)]

    values_combined = pd.concat(values)
    values_combined.drop_duplicates(subset=["sample_id", "category"], inplace=True)
    values_df = values_combined.pivot(index=["sample_id","sample_name"], columns="category",values="value").reset_index()
    return values_df
    
    
    
def main():
    args = argument_parser()
    args = argument_parser()
    username = os.getenv("GEMMA_USERNAME")
    password = os.getenv("GEMMA_PASSWORD")
    client = gemmapy.GemmaPy(auth=[username, password], path='dev')
    study_name = args.study_name
    samples_raw = client.raw.get_dataset_samples(study_name)
    dataset_info = client.get_dataset_annotations(study_name)
    
    # store assay type and platform
    # try except
    try:
        assay = dataset_info[dataset_info["class_name"]=="assay"]["term_name"].values
        assay=assay[0]
    except IndexError:
        assay = "unknown"
    samples = client.get_dataset_samples(study_name, use_processed_quantitation_type=False)
    sample_names = [x for x in samples["sample_name"]]
    sample_ids = [x.id for x in samples_raw.data]
    organisms = [x.array_design.taxon.scientific_name.lower().replace(" ", "_") for x in samples_raw.data]
    organism = organisms[0]
    
    # combine sample factor values and sample characteristics
    factor_values = [df for df in samples["sample_factor_values"]]
    sample_characteristics = [df for df in samples["sample_characteristics"]]
    # add sample id to each df in sample meta
    
    factor_meta_df = extract_meta(factor_values, sample_ids, sample_names)
    sample_meta_df = extract_meta(sample_characteristics, sample_ids, sample_names)
    
    unique_cols = [col for col in factor_meta_df.columns if col not in sample_meta_df.columns]
    # left join on sample_id
    sample_meta_df = sample_meta_df.merge(factor_meta_df[["sample_id"] + unique_cols], on=["sample_id"], how="left")
    
    
    sample_meta_df["assay"] = assay
    sample_meta_df["organism"] = organism
    

    outdir = organisms[0]
    os.makedirs(outdir, exist_ok=True)

    sample_meta_df.to_csv(f"{study_name}_sample_meta.tsv", index=False, sep="\t")
    

    
if __name__ == "__main__":
    main()
