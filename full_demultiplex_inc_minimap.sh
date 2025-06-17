#!/usr/bin/env bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <barcodes.fa> <reads.fastq.gz> <output_directory>"
    echo -e "\n<barcodes.fa> should contain the full dual barcodes with the names of the barcodes as the headers"
    echo "should be the ONT barcode + adapter + inner barcode + primer"
    echo -e "\nthis script assumes the length of the full dual barcode is ~50 bp, if the length is very different you may want to adjust the minimap2 arguments in the script \(particularly the -m argument is set to the length of the barcode minus 10\)"

    exit 1
fi
BARCODES_FA="$1"
READS_FASTQ="$2"
OUTPUT_DIR="$3"


OUTPUT_FILE="barcode_assignments.paf"


mkdir -p "$OUTPUT_DIR"
echo "Created output directory: $OUTPUT_DIR"
mkdir "$OUTPUT_DIR/intermediate_files"
# Define the log file path
LOGFILE="$OUTPUT_DIR/run.log"

# Redirect all stdout and stderr to both console and log file
exec > >(tee "$LOGFILE") 2>&1

echo -e "\nRunning minimap2..."

minimap2 "$BARCODES_FA" "$READS_FASTQ" -k5 -A1 -m32 -w1 > "$OUTPUT_DIR/intermediate_files/barcode-read_alignment.paf"

echo "minimap2 finished"

INPUT_FILE="$OUTPUT_DIR/intermediate_files/barcode-read_alignment.paf"
# Count total rows
TOTAL_ROWS=$(wc -l < "$INPUT_FILE")

# Filter and sort:
# keep reads with barcodes at the first or last 10% of the read
# keep only reads between 3 and 9 kb
# sort by number of bases matching reference barcode, descending
echo -e "\nremoving hits from minimap2 output that are from reads that are not between 3 and 9 kb, and removing hits that have barcodes in the middle"
FILTERED=$(awk -F'\t' '$2 >= 3000 && $2 <= 9000 && (($3 / $2) < 0.1 || ($3 / $2) > 0.9)' "$INPUT_FILE" | sort -k10,10nr)

# keep only the first appearance of each read ID (col 1). at this point they are sorted by number of bases matched to the barcode so this should keep the best match for each read if there are duplicates
DEDUPED=$(echo "$FILTERED" | awk -F'\t' '!seen[$1]++')

# Write to output file
echo "$DEDUPED" > "$OUTPUT_DIR/intermediate_files/$OUTPUT_FILE"


# Count kept and removed
KEPT_ROWS=$(echo "$DEDUPED" | wc -l)
REMOVED_ROWS=$((TOTAL_ROWS - KEPT_ROWS))

# Report
echo  "Filtered and sorted $INPUT_FILE -> $OUTPUT_FILE"
echo "Removed $REMOVED_ROWS rows"
echo "Remaining $KEPT_ROWS rows are written to $OUTPUT_FILE"

mkdir "$OUTPUT_DIR/intermediate_files/sorted_untrimmed_fastq_files"

# Create empty fastq files for each barcode
created_files=()
while IFS= read -r line; do
    if [[ $line == \>* ]]; then
        barcode_name=$(echo "${line#>}" | tr -d '\r')
        barcode_file="${OUTPUT_DIR}/intermediate_files/sorted_untrimmed_fastq_files/${barcode_name}.fastq"
        touch "$barcode_file"
        created_files+=("$barcode_file")
    fi
done < "$BARCODES_FA"

# Create additional fastq for unclassified reads
unclassified_file="$OUTPUT_DIR/intermediate_files/sorted_untrimmed_fastq_files/unclassified.fastq"
touch "$unclassified_file"
created_files+=("$unclassified_file")
echo -e "\nCreated empty fastq files for each barcode"

#Create an associative array from TSV (read ID -> barcode)
declare -A READ_TO_BARCODE
while IFS=$'\t' read -r read_id _ _ _ _ barcode _; do
    READ_TO_BARCODE["$read_id"]="$barcode"
done < <(cat "$OUTPUT_DIR/intermediate_files/$OUTPUT_FILE")


#Read the gzipped FASTQ and demultiplex
echo -e "\nSorting reads from $READS_FASTQ based on minimap2 output..."

zcat "$READS_FASTQ" | paste - - - - | while IFS=$'\t' read -r header seq plus qual; do
    read_id=$(echo "$header" | grep -oP 'parent_read_id=\K[^\s]+')
    barcode="${READ_TO_BARCODE[$read_id]}"
    if [ -n "$barcode" ]; then
        echo -e "@$read_id\n$seq\n+\n$qual" >> "$OUTPUT_DIR/intermediate_files/sorted_untrimmed_fastq_files/${barcode}.fastq"
    else
        echo -e "@$read_id\n$seq\n+\n$qual" >> "$OUTPUT_DIR/intermediate_files/sorted_untrimmed_fastq_files/unclassified.fastq"
    fi
done

echo "All reads have been sorted into barcode-specific fastq files."

echo -e "\nTrimming off barcodes and adapters with porechop..."
echo "IMPORTANT: for this step to work, the custom inner barcodes and primers must be added as adapters in porechop's adapters.py file"
mkdir "$OUTPUT_DIR/sorted_trimmed_fastq_files"
touch "$OUTPUT_DIR/intermediate_files/porechop.log"
for fastq in "$OUTPUT_DIR/intermediate_files/sorted_untrimmed_fastq_files"/*.fastq; do
    filename=$(basename "$fastq")
    echo "Chopping $filename"
    echo -e "\nStarting $filename" >> "$OUTPUT_DIR/intermediate_files/porechop.log"

    porechop -i "$fastq" -o "$OUTPUT_DIR/sorted_trimmed_fastq_files/trimmed_$filename" --discard_middle >>"$OUTPUT_DIR/intermediate_files/porechop.log"


    echo "Finished $filename" >>"$OUTPUT_DIR/intermediate_files/porechop.log"
done
echo -e "\n porechop finished (find detailed info in $OUTPUT_DIR/intermediate_files/porechop.log"
echo -e "\nFinal Barcode distribution:"

# Loop over each FASTQ file in the output folder and generate ASCII bars
for fastq in "$OUTPUT_DIR/sorted_trimmed_fastq_files"/*.fastq; do
    count=$(wc -l < "$fastq")
    reads=$((count / 4))
    bar=$(printf '%*s' "$((reads / 100))" '' | tr ' ' '#')  # Adjust divisor to scale
    filename=$(basename "$fastq")
    printf "%-20s %5d %s\n" "$filename"          "$reads" "$bar"
done

echo -e "\nAll the above output can also be found in $OUTPUT_DIR/run.log"
