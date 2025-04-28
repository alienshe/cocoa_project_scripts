#!/usr/bin/env bash


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <dir/of/devider/pipline/outputs> <output_file.txt>"
    exit 1
fi

PARENT_DIR="$1"
OUTPUT="$2"

> "$OUTPUT"

for dev_out in "$PARENT_DIR"/*; do
    maj_vote_file="$dev_out/devider_output/majority_vote_haplotypes.fasta"
    base_folder_name=$(basename "$dev_out")
    echo -e "\n$base_folder_name" >> "$OUTPUT"
    if [ -f "$maj_vote_file" ]; then
        grep "^>" "$maj_vote_file" >> "$OUTPUT"
    else
        echo "No file found: $maj_vote_file" >> "$OUTPUT"
    fi

done

cat "$OUTPUT"


