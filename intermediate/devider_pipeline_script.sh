#!/bin/bash
GENE_NAME="Tc3618"
INPUT_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9/amplicons/Tc3618-amplicons/porechopped/Tc3618_reads_porechopped.fastq"
PARENT_OUTPUT_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9/devider_output"
REFERENCE_FILE_NAME=""
REFERENCE_FILE_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9/reference_files"
THREADS="8"
OUTPUT_DIR="${PARENT_OUTPUT_DIR}/${GENE_NAME}_pipeline_porechopped_reads_changed_headers"
touch ~/sylvie_workspace/temp_log.txt
run_devider_pipeline \
-i ${INPUT_DIR} \
-r ${REFERENCE_FILE_DIR}/${GENE_NAME}_ref_chr9.fasta \
-o ${OUTPUT_DIR} \
-t ${THREADS} \
--overwrite \
>~/sylvie_workspace/temp_log.txt 2>&1

mv ~/sylvie_workspace/temp_log.txt ${OUTPUT_DIR}/log.txt

