#!/usr/bin/env python3
"""Rename metadata TSV columns to a standardized schema."""

import argparse
import pandas as pd

RENAME_MAP = {
    # sex
    "Biological_Sex": "sex",
    "biological sex": "sex",
    "Sex": "sex",
    "msex": "sex",
    # age
    "Age_death": "age",
    "age_death": "age",
    "Age": "age",
    "Age of death (Y)": "age",
    # pmi
    "PMI": "pmi",
    "PMI_h": "pmi",
    # region
    "organism part": "region",
    # disease
    "Disorder": "disease",
    "diagnosis": "disease",
    "Condition": "disease",
    "Pathologic_diagnosis_of_AD": "disease",
    "Schizophrenia": "disease",
    # donor_id
    "Individual_ID": "donor_id",
    "Donor": "donor_id",
    "subject": "donor_id",
    "Mayo_ID": "donor_id",
    # ethnicity
    "1000G_ancestry": "ethnicity",
    "Race": "ethnicity",
    "race": "ethnicity",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input metadata TSV")
    parser.add_argument("-o", "--output", required=True, help="Output metadata TSV")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    renames = {col: RENAME_MAP[col] for col in df.columns if col in RENAME_MAP}
    if renames:
        df.rename(columns=renames, inplace=True)
        for old, new in sorted(renames.items()):
            print(f"  Renamed: {old} -> {new}")
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
