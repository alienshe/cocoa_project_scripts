#!/bin/bash

#positional args:
INPUT_DIR=$1 #directory containing reads and ref fasta, absolute path needed
READS=$2 #relative to input directory, must be trimmed already eg with porechop
REFERENCE_FASTA=$3 #relative to input directory, subset of the relevant chromosome, no special characters in header of fasta, must be indexed with samtools (.fasta.fai file in same directory)
OUTPUT_DIR=$4 #absolute path needed
MINIMUM_ABUNDANCE=${5:-"0"} #optional, number from 0-100,(default is 0) which gets inputted to devider 
GENE_NAME=${6:-$(basename "${OUTPUT_DIR}")} #optional, use if you want intermediate files like alignments to be named something different than your output directory

#use for one gene at a time

#example: bash sylvies_devider_pipeline /home/user/upper_level_dir my_amplicon_folder reference_files/chr9.fasta /home/user/output Tc3618
LOG_FILE="${OUTPUT_DIR}/log.txt"
basename "${OUTPUT_DIR}"
#other settings
THREADS="8" #used for clair3 and devider
MODEL_NAME="r1041_e82_400bps_hac_v500"
MODEL_DIR="/home/tobias-lab/models"


#make output_dir
echo "Creating output directory"
mkdir ${OUTPUT_DIR}
echo "Starting pipeline for gene: ${GENE_NAME}" | tee ${LOG_FILE}  # This ensures you log the start of the pipeline
exec > >(tee -a ${LOG_FILE}) 2>&1  # Redirect both stdout and stderr to the log file

echo $@>>${OUTPUT_DIR}/log.txt #write all arguments to the log
echo "Model used: "${MODEL_NAME} >>${OUTPUT_DIR}/log.txt


#start making alignments
echo "Starting to constuct alignments with minimap"
mkdir ${OUTPUT_DIR}/alignments
cd ${OUTPUT_DIR}/alignments
#make initial alignment
minimap2 -ax map-ont\
 ${INPUT_DIR}/${REFERENCE_FASTA}\
 ${INPUT_DIR}/${READS}\
 > ${GENE_NAME}_alignment.sam
#convert to bam file
samtools view -bS\
 ${GENE_NAME}_alignment.sam\
 > ${GENE_NAME}_alignment.bam

samtools sort\
 -o ${GENE_NAME}_alignment_sorted.bam\
 ${GENE_NAME}_alignment.bam

samtools index\
 -o ${GENE_NAME}_alignment_sorted.bam.bai\
 ${GENE_NAME}_alignment_sorted.bam

echo "Finished making alignments"

#start making vcf with clair3

echo -e "\nStarting vcf production with clair3"
echo "alignment: ${OUTPUT_DIR}/alignments/${GENE_NAME}_alignment_sorted.bam"
echo "reference: ${INPUT_DIR}/${REFERENCE_FASTA}"
mkdir ${OUTPUT_DIR}/clair3_output
docker run \
-v ${MODEL_DIR}:/opt/models \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
hkubal/clair3:latest /opt/bin/run_clair3.sh \
--bam_fn=${OUTPUT_DIR}/alignments/${GENE_NAME}_alignment_sorted.bam \
--ref_fn=${INPUT_DIR}/${REFERENCE_FASTA} \
--threads=${THREADS} \
--platform="ont" \
--model_path="/opt/models/${MODEL_NAME}" \
--output=${OUTPUT_DIR}/clair3_output \
--include_all_ctgs \
--chunk_size=5000 \
--enable_variant_calling_at_sequence_head_and_tail \
--enable_phasing \
--var_pct_full=1.0 \
--ref_pct_full=1.0
echo "Finished clair3"

#devider

echo "Starting devider"

mkdir ${OUTPUT_DIR}/devider_output

devider \
-v ${OUTPUT_DIR}/clair3_output/merge_output.vcf.gz \
-b ${OUTPUT_DIR}/alignments/${GENE_NAME}_alignment_sorted.bam \
-r ${INPUT_DIR}/${REFERENCE_FASTA} \
-o ${OUTPUT_DIR}/devider_output \
-t ${THREADS} \
--allele-output \
--preset nanopore-r10 \
--min-abund ${MINIMUM_ABUNDANCE}

haplotag_bam \
${OUTPUT_DIR}/alignments/${GENE_NAME}_alignment_sorted.bam \
-i ${OUTPUT_DIR}/devider_output/ids.txt

echo "Finished devider"


