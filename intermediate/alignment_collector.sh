#!/usr/bin/env bash


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <dir/of/devider/pipline/outputs> <output/dir.txt>"
    exit 1
fi

PARENT_DIR="$1"
OUTPUT="$2"

mkdir "$OUTPUT"
mkdir ${OUTPUT}/whatshap
mkdir ${OUTPUT}/whatshap/individual_consensus_seq_files
mkdir ${OUTPUT}/whatshap/alignments
mkdir ${OUTPUT}/devider
parent_basename=$(basename "$PARENT_DIR")
for dev_out in "$PARENT_DIR"/*; do
    if [ ! -d "$dev_out"/alignments ]; then
        continue
    fi

    base_folder_name=$(basename "$dev_out")
    alignment_dir=$"$dev_out"/alignments
    if [ -f "$alignment_dir"/*_alignment_sorted.bam.tagged.bam ]; then
        cp -r "$alignment_dir"/"$base_folder_name"_alignment_sorted.bam.tagged.bam "$OUTPUT"/devider/"$base_folder_name"_"$parent_basename"_alignment_sorted_tagged_devider.bam
        cp -r "$alignment_dir"/"$base_folder_name"_alignment_sorted.bam.tagged.bam.bai "$OUTPUT"/devider/"$base_folder_name"_"$parent_basename"_alignment_sorted_tagged_devider.bam.bai
    else
        cp -r "$alignment_dir"/"$base_folder_name"_alignment_sorted.bam "$OUTPUT"/devider/"$base_folder_name"_"$parent_basename"_alignment_sorted.bam
        cp -r "$alignment_dir"/"$base_folder_name"_alignment_sorted.bam.bai "$OUTPUT"/devider/"$base_folder_name"_"$parent_basename"_alignment_sorted.bam.bai
    fi
    cp ${dev_out}/whatshap_output/* ${OUTPUT}/whatshap
    mv ${OUTPUT}/whatshap/*.fasta ${OUTPUT}/whatshap/individual_consensus_seq_files/
    mv ${OUTPUT}/whatshap/*.bam ${OUTPUT}/whatshap/alignments/

done

for consensus_seq in ${OUTPUT}/whatshap/individual_consensus_seq_files/*; do
    original_header=$(awk 'sub(/^>/, "")' $consensus_seq)
    new_header=$(basename $consensus_seq)
    touch ${OUTPUT}/whatshap/${original_header}_consensus_sequences.fasta
    awk -v new_header="$new_header" 'FNR==1{print ">" new_header; next}1' $consensus_seq >> ${OUTPUT}/whatshap/${original_header}_consensus_sequences.fasta
done
