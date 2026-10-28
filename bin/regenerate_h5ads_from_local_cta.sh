#!/bin/bash
set -uo pipefail

# Regenerates per-sample h5ads for one or more studies using a local
# cell_type_assignments/<study>.celltypes.tsv instead of Gemma's, without
# going through Nextflow. For studies where Gemma can't yet serve the
# correct CTA (e.g. blocked on curation), drop the correct file into
# cell_type_assignments/ first, then run this.
#
# Usage: regenerate_h5ads_from_local_cta.sh <outdir> <study_name> [study_name ...]

outdir=$1
shift
studies=("$@")

if [ -z "$outdir" ] || [ ${#studies[@]} -eq 0 ]; then
    echo "Usage: $0 <outdir> <study_name> [study_name ...]"
    exit 1
fi

gene_map=/space/grp/rschwartz/rschwartz/get_gemma_data.nf/meta/gemma_genes.tsv
py=/home/rschwartz/anaconda3/envs/scanpyenv/bin/python
script=/space/grp/rschwartz/rschwartz/get_gemma_data.nf/bin/process_query_samples.py
work_dir=$(mktemp -d)
log="$outdir/regenerate_h5ads.log"
> "$log"

for study in "${studies[@]}"; do
    for sample_dir in "$outdir/mex/$study"/*/; do
        sample_name=$(basename "$sample_dir")
        run_dir="$work_dir/${study}_${sample_name}"
        mkdir -p "$run_dir"
        cd "$run_dir" || exit 1
        echo "=== $study $sample_name ===" >> "$log"
        "$py" "$script" \
            --query_path "$sample_dir" \
            --cell_meta_path "$outdir/cell_type_assignments/${study}.celltypes.tsv" \
            --sample_meta_path "$outdir/metadata_standardized/$study/${study}_sample_meta_std.tsv" \
            --query_name "$sample_name" \
            --gene_mapping "$gene_map" \
            --study_name "$study" >> "$log" 2>&1

        if [ -f "${study}_${sample_name}.h5ad" ]; then
            mkdir -p "$outdir/h5ad/$study"
            cp "${study}_${sample_name}.h5ad" "$outdir/h5ad/$study/"
            echo "  -> updated h5ad/$study/${study}_${sample_name}.h5ad" >> "$log"
        elif [ -f "small_samples/${sample_name}.h5ad" ]; then
            mkdir -p "$outdir/small_samples/$study"
            cp "small_samples/${sample_name}.h5ad" "$outdir/small_samples/$study/"
            echo "  -> moved to small_samples/$study/${sample_name}.h5ad (fewer than 50 verified cells)" >> "$log"
            rm -f "$outdir/h5ad/$study/${study}_${sample_name}.h5ad"
        else
            echo "  -> NO OUTPUT PRODUCED (unexpected)" >> "$log"
        fi
        cd - > /dev/null
    done
done

rm -rf "$work_dir"
echo "DONE" >> "$log"
