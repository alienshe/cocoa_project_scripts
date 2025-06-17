#!/usr/bin/env bash


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <dir/of/devider/pipline/outputs> <output_file.txt>"
    exit 1
fi

PARENT_DIR="$1"
OUTPUT="$2"

> "$OUTPUT"

for dev_out in "$PARENT_DIR"/*; do
    if [ -d "$dev_out/alignments" ]; then
        #find haplotype file (maj_vote_file)
        maj_vote_file="$dev_out/devider_output/majority_vote_haplotypes.fasta"
        #print name of the folder (should be name of the gene)
        base_folder_name=$(basename "$dev_out")
        echo -e "\n$base_folder_name" >> "$OUTPUT"
        if [ -f "$maj_vote_file" ]; then
            grep "^>" "$maj_vote_file" >> "$OUTPUT"
        else
            echo "No devider output found, probably means that no variant SNPs were detected by Clair3" >> "$OUTPUT"
        fi
    else
        continue
    fi

done

cat "$OUTPUT"


