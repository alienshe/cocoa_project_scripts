#! /bin/bash
source /home/tobias-lab/miniconda3/bin/activate whatshap-env
OUTPUT_DIRECTORY="$1"
OUTPUT_ALIGNMENT="$2"
REFERENCE="$3"
PHASED_VCF="$4"
INPUT_ALIGNMENT="$5"
GENE_NAME="$6"

echo ""
echo ""
echo ""
echo "haplotagging with whatshap"
whatshap haplotag --ignore-read-groups \
	-o "$OUTPUT_DIRECTORY"/"$OUTPUT_ALIGNMENT"  \
	--reference $REFERENCE \
	$PHASED_VCF \
	$INPUT_ALIGNMENT

source /home/tobias-lab/miniconda3/bin/activate base
echo ""
echo ""
echo ""
echo "splitting reads based on whatshap output"
mkdir ${OUTPUT_DIRECTORY}/split_alignments
samtools split \
	-d HP \
	-u ${OUTPUT_DIRECTORY}/split_alignments/${GENE_NAME}_haplotype_unassigned_alignment.bam \
	-f "${OUTPUT_DIRECTORY}/split_alignments/${GENE_NAME}_haplotype_%#_alignment.%." \
	"$OUTPUT_DIRECTORY"/"$OUTPUT_ALIGNMENT"


echo ""
echo ""
echo ""
echo "creating consensus sequence for each haplotype"
for alignment in ${OUTPUT_DIRECTORY}/split_alignments/*; do
    
	base_name=$(basename "$alignment")
	stem=$(grep -Po '.*_haplotype_.*_' <<< "$base_name")
	samtools consensus -o "${OUTPUT_DIRECTORY}/${stem}consensus.fasta" $alignment 
done
