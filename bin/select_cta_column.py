#!/usr/bin/env python3
"""Convert a `gemma-cli-staging getSingleCellMetadata -allQts` TSV into the
celltypes.tsv schema the rest of the pipeline expects: sample_id, cell_id,
cell_type, cell_type_uri.

Background: `getSingleCellMetadata -allQts` returns every available cell type
assignment (CTA) / quantitation type for a dataset as its own column, instead
of guessing one protocol name up front like the curl-based REST call does.
Those CTA columns are appended after the fixed cellId/Bioassay/technical
columns and whatever sample-characteristic columns the dataset happens to
have, so which columns are "real" CTAs isn't knowable from the header alone
in general. This script uses a simple default for picking one:

  1. Prefer a column whose name looks like an author-submitted CTA (contains
     "author", case-insensitively) -- this matches the common
     "author-submitted" protocol name used elsewhere in this pipeline.
  2. Otherwise, fall back to the last column in the file, since Gemma
     appends CTA/quantitation-type columns after everything else.

This is a heuristic, not a guarantee -- it's logged clearly so the choice can
be sanity-checked per dataset.

Expects the input to have been produced with:
    gemma-cli-staging getSingleCellMetadata -e <study> -allQts \\
        -useBioAssayIds -useRawColumnNames -o <input>
`-useBioAssayIds` makes `cellId` take the form "<bioassayId>_<barcode>",
where <bioassayId> matches the numeric prefix used for MEX sample directory
names elsewhere in this pipeline, and <barcode> matches the raw 10x barcode
(e.g. "AAACCTGAGATAGCAT-1") used as the AnnData obs index. `-useRawColumnNames`
keeps CTA protocol names intact (e.g. "author-submitted" instead of
"author.submitted") and preserves "-" in barcodes instead of converting to ".".
"""

import argparse
import re

import pandas as pd

# Columns produced by `getSingleCellMetadata -allQts -useBioAssayIds
# -useRawColumnNames` that are never cell-type-assignment columns.
KNOWN_NON_CTA_COLUMNS = {
    "cellId",
    "Bioassay",
    "Sample",
    "Assays",
    "number of cells",
    "number of cells by design elements",
    "number of design elements",
    "sequence read count",
}


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Raw TSV from `gemma-cli-staging getSingleCellMetadata -allQts`")
    parser.add_argument("--output", required=True, help="Output celltypes TSV path (sample_id/cell_id/cell_type/cell_type_uri)")
    parser.add_argument("--study_name", required=True, help="Study name, used only for log messages")
    return parser.parse_args()


def select_cta_column(columns, study_name):
    candidates = [c for c in columns if c not in KNOWN_NON_CTA_COLUMNS]
    if not candidates:
        raise ValueError(f"[{study_name}] No candidate cell type assignment columns found among: {list(columns)}")

    author_like = [c for c in candidates if re.search(r"author", c, re.IGNORECASE)]
    if author_like:
        chosen = author_like[0]
        if len(author_like) > 1:
            print(f"[{study_name}] Multiple author-like CTA columns found {author_like}; using '{chosen}'.")
    else:
        chosen = candidates[-1]

    print(f"[{study_name}] Candidate CTA/quantitation-type columns: {candidates}")
    print(f"[{study_name}] Selected CTA protocol column: '{chosen}'")
    return chosen


def main():
    args = parse_arguments()
    df = pd.read_csv(args.input, sep="\t")

    missing = {"cellId", "Bioassay"} - set(df.columns)
    if missing:
        raise ValueError(f"[{args.study_name}] Input {args.input} is missing expected column(s) {missing}: {list(df.columns)}")

    cta_column = select_cta_column(df.columns, args.study_name)

    # With -useBioAssayIds, cellId is "<bioassayId>_<barcode>"; <bioassayId>
    # matches the numeric prefix of MEX sample directory names, and <barcode>
    # matches the raw AnnData obs index used downstream.
    split_id = df["cellId"].astype(str).str.split("_", n=1, expand=True)
    if split_id.shape[1] < 2 or split_id[1].isna().any():
        raise ValueError(
            f"[{args.study_name}] Could not split 'cellId' into <bioassayId>_<barcode> for all rows; "
            f"was the input generated with -useBioAssayIds? Example values: {df['cellId'].head(3).tolist()}"
        )

    out = pd.DataFrame({
        "sample_id": split_id[0],
        "cell_id": split_id[1],
        "cell_type": df[cta_column],
        "cell_type_uri": "",
    })
    out.to_csv(args.output, sep="\t", index=False)
    print(f"[{args.study_name}] Wrote {len(out)} rows to {args.output}")


if __name__ == "__main__":
    main()
