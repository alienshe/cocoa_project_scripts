#!/bin/bash

INPUT_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="/home/tobias-lab/sylvie_workspace/NLR_HD_amplicons/chr9/clair3_output"      # e.g. /home/user1/output (absolute path needed)
THREADS="8"            # e.g. 8
MODEL_NAME="r941_prom_hac_g360+g422"
#REFERENCE_FILE_NAME=${GENE_NAME}_ref_chr9
GENE_NAME="Tc3425"
source activate base
conda activate clair3
docker run -it \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
hkubal/clair3:latest /opt/bin/run_clair3.sh \
--bam_fn=${INPUT_DIR}/alignments/${GENE_NAME}_alignment_sorted.bam \ #input bam file
--ref_fn=${INPUT_DIR}/reference_files/${GENE_NAME}_ref_chr9.fasta \ #reference file
--threads=${THREADS} --platform="ont" \ 
--model_path="/opt/models/${MODEL_NAME}" \
--output=${OUTPUT_DIR}/${GENE_NAME} \
--include_all_ctgs  #recommended for non-human genomes
