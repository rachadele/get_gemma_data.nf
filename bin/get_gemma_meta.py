#!/bin/python
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import argparse
import gemmapy

def argument_parser():
    parser = argparse.ArgumentParser(description="Preprocess data from GEMMA")
    parser.add_argument("--study_name", type=str, help="Name of the study", default="GSE152715.1")
    parser.add_argument('--gemma_username', type=str, required=True)
    parser.add_argument('--gemma_password', type=str, required=True)
    return parser.parse_args()

def main():
    args = argument_parser()

    client = gemmapy.GemmaPy(auth=[args.gemma_username, args.gemma_password], path='staging')
    study_name = args.study_name

    # use_processed_quantitation_type=False is the raw per-BioAssay set this
    # pipeline needs (matches MEX/h5ad filenames, CTA files); each object
    # already carries its own BioMaterial (`.sample`) with characteristics
    # nested inside, so id/name/characteristics all come from one place with
    # no risk of misalignment. (=True instead explodes each sample into one
    # BioAssay per cell type, unrelated to what's needed here.)
    samples_raw = client.raw.get_dataset_samples(study_name, use_processed_quantitation_type=False)

    # Some datasets (e.g. SZBDMulti-Seq, GSE254569) come back with a null
    # array_design.taxon on every BioAssay; fall back to the dataset's own taxon
    # (looked up lazily, only if needed) and fail loudly if that is missing too.
    dataset_taxon = None

    def organism_of(s):
        nonlocal dataset_taxon
        if s.array_design is not None and s.array_design.taxon is not None:
            return s.array_design.taxon.scientific_name.lower().replace(" ", "_")
        if dataset_taxon is None:
            found = client.raw.get_datasets(filter=f"shortName = {study_name}").data
            if len(found) != 1 or found[0].taxon is None:
                raise ValueError(f"{study_name}: BioAssay taxon is null and no dataset-level taxon found")
            dataset_taxon = found[0].taxon.scientific_name.lower().replace(" ", "_")
        return dataset_taxon

    rows = []
    for s in samples_raw.data:
        row = {
            "sample_id": s.id,             # BioAssay ID (join key)
            "biomaterial_id": s.sample.id,  # BioMaterial ID (reference only)
            "sample_name": s.name,
            "organism": organism_of(s),
        }
        for c in s.sample.characteristics:
            if c.category is not None:
                row[c.category] = c.value
        rows.append(row)

    sample_meta_df = pd.DataFrame(rows)
    sample_meta_df.to_csv(f"{study_name}_sample_meta.tsv", index=False, sep='\t')


if __name__ == "__main__":
    main()
