#!/bin/bash

#positional args:
INPUT_VCF=$1
REFERENCE_FASTA=$2
ALIGNMENT=$3
OUTPUT_DIR=$4
MINIMUM_ABUNDANCE=${5:-"0"} #optional, number from 0-100,(default is 0) which gets inputted to devider
THREADS=${6:-"8"}


DEVIDER_PATH="/home/tobias-lab/miniconda3/bin/devider"

mkdir ${OUTPUT_DIR}

$DEVIDER_PATH \
-v ${INPUT_VCF} \
-b ${ALIGNMENT} \
-r ${REFERENCE_FASTA} \
-o ${OUTPUT_DIR} \
-t ${THREADS} \
--allele-output \
--preset nanopore-r10 \
--min-abund ${MINIMUM_ABUNDANCE}

haplotag_bam \
$ALIGNMENT \
-i ${OUTPUT_DIR}/ids.txt

echo "Finished devider"


