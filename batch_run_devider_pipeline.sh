#!/usr/bin/env bash

INPUT_TSV="$1"
INPUT_DIR="$2"
OUTPUT_DIR="$3"
QUERY_DIR="$4"
REF_DIR="$5"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Failed to cd into $OUTPUT_DIR"; exit 1; }

while IFS=$'\t' read -r query reference gene_name; do
    query=$(echo "$query" | tr -d '\r')
    reference=$(echo "$reference" | tr -d '\r')
    gene_name=$(echo "$gene_name" | tr -d '\r')
    echo "Processing:"
    echo "Query: $query"
    echo "Reference: $reference"
    echo "Gene: $gene_name"

    # Example command using query and reference
    ~/sylvie_workspace/scripts/sylvies_devider_pipeline.sh "$INPUT_DIR" "$QUERY_DIR/$query" "$REF_DIR/$reference" "$OUTPUT_DIR/$gene_name" "$gene_name"

done < "$INPUT_TSV"

echo -e "\n******SUMMARY******"
~/sylvie_workspace/scripts/devider_output_summarizer.sh "$OUTPUT_DIR" "$OUTPUT_DIR/haplotype_summary.txt"
