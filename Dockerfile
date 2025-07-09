FROM mambaorg/micromamba:1.5

USER root

# Install system deps and bioinformatics tools in one layer
RUN apt-get update && apt-get install -y git wget && rm -rf /var/lib/apt/lists/* && \
    micromamba install -y -n base -c bioconda -c conda-forge \
        porechop minimap2 samtools tabix lofreq devider \
        python pip numpy pandas biopython matplotlib seaborn && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/bin:$PATH"

# Copy files and setup in one layer
COPY . /app/
WORKDIR /app
RUN chmod +x *.sh

# Create entrypoint script
RUN cat > /usr/local/bin/entrypoint.sh << 'EOF'
#!/bin/bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate base

if [ $# -eq 0 ]; then
    cat << 'HELP'
Cocoa Project Genomics Pipeline

Command Menu:
  list_demo_files                    - List available demo files
  demo_demultiplex <output_dir>      - Run demultiplexing demo
  demo_devider <gene> <output_dir>   - Run single gene analysis demo
  demo_batch_devider <output_dir>    - Run batch processing demo
  demultiplex <args>                 - Run demultiplexing on your data
  devider <args>                     - Run single gene analysis
  batch_devider <args>               - Run batch processing

Direct Usage:
  ./script.sh <args>                 - Run any script directly
  bash                               - Interactive shell

Tools available: porechop, minimap2, samtools, devider
HELP
    command -v devider >/dev/null && echo "✅ devider ready" || echo "⚠️  devider not found"
    exit 0
fi

case "$1" in
    list_demo_files) find /app/demo_files -type f | head -20 ;;
    demo_demultiplex) 
        [ -z "$2" ] && { echo "Usage: demo_demultiplex <output_dir>"; exit 1; }
        mkdir -p "$2"
        ./full_demultiplex_inc_minimap.sh \
            demo_files/inputs/reference_files/dual_barcodes_with_primers.fa \
            demo_files/inputs/pangkep_raw_reads/barcode01_raw_reads.fastq.gz "$2" ;;
    demo_devider)
        [ -z "$2" ] || [ -z "$3" ] && { 
            echo "Usage: demo_devider <gene_name> <output_dir>"
            echo "Available genes: Tc1318, Tc1320, Tc2391, Tc3078, Tc3425, Tc3618"
            exit 1; }
        mkdir -p "$3"
        ./sylvies_devider_pipeline.sh /app/demo_files \
            "demo_files/inputs/pangkep_run_bc01_demultiplexed/trimmed_$2.fastq" \
            "demo_files/inputs/reference_files/${2}_ref_chr*.fasta" "$3" 5 "$2" ;;
    demo_batch_devider)
        [ -z "$2" ] && { echo "Usage: demo_batch_devider <output_dir>"; exit 1; }
        mkdir -p "$2"
        ./batch_run_devider_pipeline.sh \
            demo_files/inputs/reference_files/NLRs.tsv /app/demo_files "$2" \
            demo_files/inputs/pangkep_run_bc01_demultiplexed \
            demo_files/inputs/reference_files 5 ;;
    demultiplex) shift; ./full_demultiplex_inc_minimap.sh "$@" ;;
    devider) shift; ./sylvies_devider_pipeline.sh "$@" ;;
    batch_devider) shift; ./batch_run_devider_pipeline.sh "$@" ;;
    bash|sh|/bin/bash|/bin/sh) exec /bin/bash ;;
    *) exec /bin/bash "$@" ;;
esac
EOF

RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]