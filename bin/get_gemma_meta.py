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

    # Each BioAssay object from this one call already carries its own
    # BioMaterial (`.sample`) and that BioMaterial's characteristics nested
    # inside it -- so sample_id/name/characteristics all come from the same
    # object, with no risk of misalignment. A previous version of this
    # script instead fetched names/characteristics from a *second*,
    # separately-ordered API call (client.get_dataset_samples(...)) and
    # zipped the two together by list position; that second call's
    # "sample_ID" turned out to be a BioMaterial ID (a different ID space
    # than BioAssay, zero overlap), so the positional zip silently attached
    # the wrong sample's name/characteristics to a given BioAssay id in
    # several studies (confirmed on CMC: 100/100 samples affected).
    samples_raw = client.raw.get_dataset_samples(study_name)

    rows = []
    for s in samples_raw.data:
        row = {
            "sample_id": s.id,
            "sample_name": s.name,
            "organism": s.array_design.taxon.scientific_name.lower().replace(" ", "_"),
        }
        for c in s.sample.characteristics:
            row[c.category] = c.value
        rows.append(row)

    sample_meta_df = pd.DataFrame(rows)
    sample_meta_df.to_csv(f"{study_name}_sample_meta.tsv", index=False, sep='\t')


if __name__ == "__main__":
    main()
