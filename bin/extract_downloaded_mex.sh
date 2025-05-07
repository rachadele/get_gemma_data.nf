#!/bin/bash

work_dir=$1
new_dir=$2

if [ -z "$work_dir" ]; then
    echo "Usage: $0 <work_dir>"
    exit 1
fi

# Create an array from unique full paths (third-to-last directory)
readarray -t mex_dirs < <(
    find "$work_dir" -name "*mtx.gz" | while read -r filepath; do
        fullpath=$(realpath "$filepath")
        dir_three_up=$(dirname "$(dirname "$fullpath")")
        echo "$dir_three_up"
    done | sort -u
)

# move to a new directory
mkdir -p "$new_dir"
for dir in "${mex_dirs[@]}"; do
    cp -r $dir $new_dir
done

