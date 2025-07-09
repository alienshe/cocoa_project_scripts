#!/usr/bin/env bash
if [ "$#" -ne 6 ]; then
    echo "Error: incorrect number of arguments"
    echo "Usage: ./batch_run_devider_pipeline.sh <list.tsv> <input_dir> <output_dir> <query_dir> <reference_dir> <minimum_abundance>"
    echo "note that input_dir and output_dir must be absolute paths (starting from home), and query_dir and reference_dir should be relative paths from input_dir"
    echo "see inside script for more information on arguments"
    exit 1
fi
INPUT_TSV="$1"
#^ tab sep file with the following columns: name of input file in query dir, name of reference file in ref dir, and the name of the gene (just so downstream file names don't get too messy
INPUT_DIR="$2"
#^ Directory containing both the query and ref dirs, must be a non-relative path eg: /home/tobias-lab/sylvie_workspace
OUTPUT_DIR="$3"
#^output directory, must be a non-relative path eg: /home/tobias-lab/devider_outputs
QUERY_DIR="$4"
#^ relative path from INPUT_DIR. Contains all input fastq files eg pangkep_1/demultiplexed_fastqs
REF_DIR="$5"
#^ relative path from INPUT_DIR. Contains all the reference files eg reference_files
MIN_AB="$6"
#^ number from  0-100  minimum abundance for a haplotype - input to devider

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Failed to cd into $OUTPUT_DIR"; exit 1; }

while IFS=$'\t' read -r query reference gene_name || [ -n "$query" ]; do
    query=$(echo "$query" | tr -d '\r')
    reference=$(echo "$reference" | tr -d '\r')
    gene_name=$(echo "$gene_name" | tr -d '\r')
    echo "Processing:"
    echo "Query: $query"
    echo "Reference: $reference"
    echo "Gene: $gene_name"

    # Example command using query and reference
    ~/scripts/sylvies_scripts/sylvies_devider_pipeline.sh "$INPUT_DIR" "$QUERY_DIR/$query" "$REF_DIR/$reference" "$OUTPUT_DIR/$gene_name" "$MIN_AB" "$gene_name"
done < "$INPUT_TSV"

echo -e "copying alignments into one folder..."
~/scripts/sylvies_scripts/intermediate/alignment_collector.sh "$OUTPUT_DIR" "$OUTPUT_DIR"/collected_alignments


echo -e "\n****** FINAL SUMMARY ******"
~/scripts/sylvies_scripts/intermediate/devider_output_summarizer.sh "$OUTPUT_DIR" "$OUTPUT_DIR/haplotype_summary.txt"
echo -e "\nAlignments collected in $OUTPUT_DIR/collected_alignments (where devider ran, alignments are tagged with haplotype 'HP')"
