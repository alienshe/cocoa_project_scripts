#!/usr/bin/env bash


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <dir/of/devider/pipline/outputs> <output/dir.txt>"
    exit 1
fi

PARENT_DIR="$1"
OUTPUT="$2"

mkdir "$OUTPUT"
parent_basename=$(basename "$PARENT_DIR")
for dev_out in "$PARENT_DIR"/*; do
    if [ ! -d "$dev_out"/alignments ]; then
        continue
    fi

    base_folder_name=$(basename "$dev_out")
    alignment_dir=$"$dev_out"/alignments
    if [ -f "$alignment_dir"/*_alignment_sorted.bam.tagged.bam ]; then
        cp -r "$alignment_dir"/"$base_folder_name"_alignment_sorted.bam.tagged.bam "$OUTPUT"/"$base_folder_name"_"$parent_basename"_alignment_sorted_tagged.bam
        cp -r "$alignment_dir"/"$base_folder_name"_alignment_sorted.bam.tagged.bam.bai "$OUTPUT"/"$base_folder_name"_"$parent_basename"_alignment_sorted_tagged.bam.bai
    else
        cp -r "$alignment_dir"/"$base_folder_name"_alignment_sorted.bam "$OUTPUT"/"$base_folder_name"_"$parent_basename"_alignment_sorted.bam
        cp -r "$alignment_dir"/"$base_folder_name"_alignment_sorted.bam.bai "$OUTPUT"/"$base_folder_name"_"$parent_basename"_alignment_sorted.bam.bai
    fi
done
