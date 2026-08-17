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

    # Need both calls: `samples` has the characteristics/metadata we want, but
    # its "sample_ID" column is a BioMaterial ID -- a different Gemma
    # table/ID-space than the BioAssay ID this pipeline joins on elsewhere
    # (MEX/h5ad filenames, CTA files); the two ID spaces have zero overlap.
    # `samples_raw` is only used to get that real BioAssay id, matched back
    # onto `samples` by sample name (BioAssay.name == its nested
    # BioMaterial.name) rather than list position, since the two calls sort
    # samples differently.
    samples = client.get_dataset_samples(study_name, use_processed_quantitation_type=False)
    samples_raw = client.raw.get_dataset_samples(study_name, use_processed_quantitation_type=False)

    bioassay_by_name = {
        s.name: {
            "bioassay_id": s.id,
            "organism": s.array_design.taxon.scientific_name.lower().replace(" ", "_"),
        }
        for s in samples_raw.data
    }

    rows = []
    for name, biomaterial_id, characteristics in zip(
        samples["sample_name"],
        samples["sample_ID"],
        samples["sample_characteristics"],
    ):
        row = {
            "sample_id": bioassay_by_name[name]["bioassay_id"],  # BioAssay ID (join key)
            "biomaterial_id": biomaterial_id,                    # BioMaterial ID (reference only)
            "sample_name": name,
            "organism": bioassay_by_name[name]["organism"],
        }
        for _, c in characteristics.iterrows():
            if pd.notna(c["category"]):
                row[c["category"]] = c["value"]
        rows.append(row)

    sample_meta_df = pd.DataFrame(rows)
    sample_meta_df.to_csv(f"{study_name}_sample_meta.tsv", index=False, sep='\t')


if __name__ == "__main__":
    main()
