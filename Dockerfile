FROM mambaorg/micromamba:1.5.10
USER root

# Install system dependencies
RUN apt-get update && apt-get install -y git wget && \
    rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools with verbose output
RUN micromamba clean --all --yes && \
    micromamba install -y -n base -c conda-forge -c bioconda --verbose \
        porechop=0.2.4 minimap2=2.24 samtools=1.15.1 tabix=1.11 lofreq=2.1.5 pigz=2.6 \
        devider clair3=0.1.12 \
        python=3.9 pip numpy pandas biopython matplotlib seaborn && \
    micromamba clean --all --yes

# Set environment path
ENV PATH="/opt/conda/bin:$PATH"

# Download Clair3 model
RUN mkdir -p /opt/models && \
    cd /opt/models && \
    wget -q https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v500.tar.gz && \
    tar -xzf r1041_e82_400bps_hac_v500.tar.gz && \
    rm r1041_e82_400bps_hac_v500.tar.gz

# Copy project files and entrypoint script
COPY . /app/
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
WORKDIR /app

# Set permissions
RUN chmod +x /app/*.sh /usr/local/bin/entrypoint.sh

# Set entrypoint
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
