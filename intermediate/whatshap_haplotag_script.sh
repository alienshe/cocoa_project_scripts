#! /bin/bash
source /home/tobias-lab/miniconda3/bin/activate whatshap-env
OUTPUT_ALIGNMENT="$1"
REFERENCE="$2"
PHASED_VCF="$3"
INPUT_ALIGNMENT="$4"

whatshap haplotag --ignore-read-groups \
	-o $OUTPUT_ALIGNMENT  \
	--reference $REFERENCE \
	$PHASED_VCF \
	$INPUT_ALIGNMENT
