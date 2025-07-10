FROM mambaorg/micromamba:1.5.8

USER root

# Install system deps and bioinformatics tools including clair3
RUN apt-get update && apt-get install -y git wget && rm -rf /var/lib/apt/lists/* && \
    micromamba install -y -n base -c bioconda -c conda-forge \
        porechop minimap2 samtools tabix lofreq devider clair3 \
        python pip numpy pandas biopython matplotlib seaborn && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/bin:$PATH"

# Download clair3 model
RUN mkdir -p /opt/models && \
    cd /opt/models && \
    wget -q https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v500.tar.gz && \
    tar -xzf r1041_e82_400bps_hac_v500.tar.gz && \
    rm r1041_e82_400bps_hac_v500.tar.gz

COPY . /app/
WORKDIR /app
RUN chmod +x *.sh

# Create entrypoint script
RUN echo '#!/bin/bash' > /usr/local/bin/entrypoint.sh && \
    echo 'eval "$(micromamba shell hook --shell bash)"' >> /usr/local/bin/entrypoint.sh && \
    echo 'micromamba activate base' >> /usr/local/bin/entrypoint.sh && \
    echo '' >> /usr/local/bin/entrypoint.sh && \
    echo 'if [ $# -eq 0 ]; then' >> /usr/local/bin/entrypoint.sh && \
    echo '    printf "%s\n" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "Cocoa Project Genomics Pipeline" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "Command Menu:" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "  list_demo_files                    - List available demo files" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "  demo_demultiplex <output_dir>      - Run demultiplexing demo" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "  demo_devider <gene> <output_dir>   - Run single gene analysis demo" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "  demo_batch_devider <output_dir>    - Run batch processing demo" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "  demultiplex <args>                 - Run demultiplexing on your data" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "  devider <args>                     - Run single gene analysis" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "  batch_devider <args>               - Run batch processing" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "Direct Usage:" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "  ./script.sh <args>                 - Run any script directly" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "  bash                               - Interactive shell" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "" \' >> /usr/local/bin/entrypoint.sh && \
    echo '        "Tools available: porechop, minimap2, samtools, devider, clair3"' >> /usr/local/bin/entrypoint.sh && \
    echo '    command -v devider >/dev/null && echo "✅ devider ready" || echo "⚠️  devider not found"' >> /usr/local/bin/entrypoint.sh && \
    echo '    command -v run_clair3.sh >/dev/null && echo "✅ clair3 ready" || echo "⚠️  clair3 not found"' >> /usr/local/bin/entrypoint.sh && \
    echo '    [ -d "/opt/models/r1041_e82_400bps_hac_v500" ] && echo "✅ clair3 model ready" || echo "⚠️  clair3 model not found"' >> /usr/local/bin/entrypoint.sh && \
    echo '    exit 0' >> /usr/local/bin/entrypoint.sh && \
    echo 'fi' >> /usr/local/bin/entrypoint.sh && \
    echo '' >> /usr/local/bin/entrypoint.sh && \
    echo 'case "$1" in' >> /usr/local/bin/entrypoint.sh && \
    echo '    list_demo_files) find /app/demo_files -type f | head -20 ;;' >> /usr/local/bin/entrypoint.sh && \
    echo '    demo_demultiplex)' >> /usr/local/bin/entrypoint.sh && \
    echo '        [ -z "$2" ] && { echo "Usage: demo_demultiplex <output_dir>"; exit 1; }' >> /usr/local/bin/entrypoint.sh && \
    echo '        mkdir -p "$2"' >> /usr/local/bin/entrypoint.sh && \
    echo '        ./full_demultiplex_inc_minimap.sh \' >> /usr/local/bin/entrypoint.sh && \
    echo '            demo_files/inputs/reference_files/dual_barcodes_with_primers.fa \' >> /usr/local/bin/entrypoint.sh && \
    echo '            demo_files/inputs/pangkep_raw_reads/barcode01_raw_reads.fastq.gz "$2" ;;' >> /usr/local/bin/entrypoint.sh && \
    echo '    demo_devider)' >> /usr/local/bin/entrypoint.sh && \
    echo '        [ -z "$2" ] || [ -z "$3" ] && {' >> /usr/local/bin/entrypoint.sh && \
    echo '            echo "Usage: demo_devider <gene_name> <output_dir>"' >> /usr/local/bin/entrypoint.sh && \
    echo '            echo "Available genes: Tc1318, Tc1320, Tc2391, Tc3078, Tc3425, Tc3618"' >> /usr/local/bin/entrypoint.sh && \
    echo '            exit 1; }' >> /usr/local/bin/entrypoint.sh && \
    echo '        mkdir -p "$3"' >> /usr/local/bin/entrypoint.sh && \
    echo '        ./sylvies_devider_pipeline.sh /app \' >> /usr/local/bin/entrypoint.sh && \
    echo '            "demo_files/inputs/pangkep_run_bc01_demultiplexed/trimmed_$2.fastq" \' >> /usr/local/bin/entrypoint.sh && \
    echo '            "demo_files/inputs/reference_files/${2}_ref_chr*.fasta" "$3" 5 "$2" ;;' >> /usr/local/bin/entrypoint.sh && \
    echo '    demo_batch_devider)' >> /usr/local/bin/entrypoint.sh && \
    echo '        [ -z "$2" ] && { echo "Usage: demo_batch_devider <output_dir>"; exit 1; }' >> /usr/local/bin/entrypoint.sh && \
    echo '        mkdir -p "$2"' >> /usr/local/bin/entrypoint.sh && \
    echo '        ./batch_run_devider_pipeline.sh \' >> /usr/local/bin/entrypoint.sh && \
    echo '            demo_files/inputs/reference_files/NLRs.tsv /app/demo_files "$2" \' >> /usr/local/bin/entrypoint.sh && \
    echo '            demo_files/inputs/pangkep_run_bc01_demultiplexed \' >> /usr/local/bin/entrypoint.sh && \
    echo '            demo_files/inputs/reference_files 5 ;;' >> /usr/local/bin/entrypoint.sh && \
    echo '    demultiplex) shift; ./full_demultiplex_inc_minimap.sh "$@" ;;' >> /usr/local/bin/entrypoint.sh && \
    echo '    devider) shift; ./sylvies_devider_pipeline.sh "$@" ;;' >> /usr/local/bin/entrypoint.sh && \
    echo '    batch_devider) shift; ./batch_run_devider_pipeline.sh "$@" ;;' >> /usr/local/bin/entrypoint.sh && \
    echo '    bash|sh|/bin/bash|/bin/sh) exec /bin/bash ;;' >> /usr/local/bin/entrypoint.sh && \
    echo '    *) exec /bin/bash "$@" ;;' >> /usr/local/bin/entrypoint.sh && \
    echo 'esac' >> /usr/local/bin/entrypoint.sh

RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]