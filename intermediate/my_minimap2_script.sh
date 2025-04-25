#!/bin/bash

GENE_NAME="Tc3425"
INPUT_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9/amplicons/${GENE_NAME}-amplicons/porechopped/${GENE_NAME}_reads_porechopped.fastq"
OUTPUT_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9/alignments/"
REFERENCE_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9/reference_files"

cd $OUTPUT_DIR

#make sam file
minimap2 -ax map-ont\
 ${REFERENCE_DIR}/${GENE_NAME}_ref_chr9.fasta\
 ${INPUT_DIR}\
 > ${GENE_NAME}_alignment.sam

#convert to bam file

samtools view -bS\
 ${GENE_NAME}_alignment.sam\
 > ${GENE_NAME}_alignment.bam

#sort
samtools sort\
 -o ${GENE_NAME}_alignment_sorted.bam\
 ${GENE_NAME}_alignment.bam

#index

samtools index\
 -o ${GENE_NAME}_alignment_sorted.bam.bai\
 ${GENE_NAME}_alignment_sorted.bam


