#!/usr/bin/env python3
"""Standardize metadata column names across h5ad files."""

import argparse
import sys
from pathlib import Path

import anndata

# Every known variant -> standardized name
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


def process_file(input_path: Path, output_path: Path) -> dict[str, str]:
    """Read an h5ad file, rename obs columns, write to output_path.

    Returns a dict of old_name -> new_name for columns that were renamed.
    """
    adata = anndata.read_h5ad(input_path)

    renames = {col: RENAME_MAP[col] for col in adata.obs.columns if col in RENAME_MAP}
    if renames:
        adata.obs.rename(columns=renames, inplace=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    return renames


def main():
    parser = argparse.ArgumentParser(
        description="Standardize metadata column names in h5ad files."
    )
    parser.add_argument(
        "input_dir",
        nargs="?",
        default="all_homo_sapiens_samples/h5ad",
        help="Input directory containing study subdirectories with h5ad files "
        "(default: all_homo_sapiens_samples/h5ad)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="all_homo_sapiens_samples/h5ad_standardized",
        help="Output directory (default: all_homo_sapiens_samples/h5ad_standardized)",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.is_dir():
        print(f"Error: input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    study_dirs = sorted(p for p in input_dir.iterdir() if p.is_dir())
    if not study_dirs:
        print(f"No study directories found in {input_dir}", file=sys.stderr)
        sys.exit(1)

    for study_dir in study_dirs:
        h5ad_files = sorted(study_dir.glob("*.h5ad"))
        if not h5ad_files:
            continue

        study_name = study_dir.name
        study_renames: dict[str, str] = {}

        for h5ad_file in h5ad_files:
            rel = h5ad_file.relative_to(input_dir)
            out_path = output_dir / rel
            print(f"  Processing {rel} ...")
            renames = process_file(h5ad_file, out_path)
            study_renames.update(renames)

        if study_renames:
            renamed_str = ", ".join(
                f"{old} -> {new}" for old, new in sorted(study_renames.items())
            )
            print(f"  [{study_name}] Renamed: {renamed_str}")
        else:
            print(f"  [{study_name}] No columns renamed")
        print()


if __name__ == "__main__":
    main()
