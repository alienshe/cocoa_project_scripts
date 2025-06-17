#!/bin/bash
INPUT_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9"
OUTPUT_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9/devider_output"
GENE_NAME="Tc3425"
THREADS="8"

devider \
-v ${INPUT_DIR}/clair3_output/${GENE_NAME}/full_alignment.vcf.gz \
-b ${INPUT_DIR}/alignments/${GENE_NAME}_alignment_sorted.bam \
-r ${INPUT_DIR}/reference_files/${GENE_NAME}_ref_chr9.fasta \
-o ${OUTPUT_DIR}/${GENE_NAME} \
-t ${THREADS}

