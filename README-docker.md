## Docker 

### Prerequisites

- Docker installed and engine running

### Pull the Docker Image

```bash  
# Pull the latest stable version  
docker pull ghcr.io/alienshe/cocoa_project_scripts:stable  

# Or pull a specific version  
docker pull ghcr.io/alienshe/cocoa_project_scripts:v1.0.0
```


### Or clone the repository and build the Docker image

```bash
git clone https://github.com/alienshe/cocoa_project_scripts.git
cd cocoa_project_scripts

# OPTIONAL - COPY THE DEMO FILES to the demo_files directory if they're needed in the image

# Default build to host architecture    
docker build -t cocoa-project .

# For ARM64 systems (Apple Silicon Mac):
docker build --platform=linux/arm64 -t cocoa-project .

# Run on Apple Silicon amd64 emulation mode:
docker build --platform=linux/amd64 -t cocoa-project .
```
## Quick Start

### Interactive Menu

```bash
# Run the container interactive shell
docker run -it cocoa-project bash

# Run the container to see the interactive menu
docker run -it cocoa-project

# Or with volume mounting for persistent data
docker run -it -v $(pwd)/data:/app/data cocoa-project
```

The interactive menu provides:

Demo commands: Try the pipeline with included test data
- Production commands: Run on your own data
- Tool access: Direct access to all bioinformatics tools
- Status checks: Verify all tools and models are properly installed

### Demo Analysis

```bash
# List available demo files
docker run -it cocoa-project list_demo_files

# Run demultiplexing demo
docker run -it -v $(pwd)/output:/app/output cocoa-project demo_demultiplex /app/output/demux_demo

# Run single gene analysis demo
docker run -it -v $(pwd)/output:/app/output cocoa-project demo_devider Tc1318 /app/output/gene_demo

# Run batch processing demo
docker run -it -v $(pwd)/output:/app/output cocoa-project demo_batch_devider /app/output/batch_demo
```

## Usage

### Demultiplexing

```bash

docker run -it \
    -v /path/to/your/data:/app/data \
    -v /path/to/output:/app/output \
    cocoa-project \
    demultiplex \
    /app/data/barcodes.fa \
    /app/data/raw_reads.fastq.gz \
    /app/output/demultiplexed
```

### Single Gene Analysis

```bash
docker run -it \
    -v /path/to/your/data:/app/data \
    -v /path/to/output:/app/output \
    cocoa-project \
    devider \
    /app/data \
    /app/data/trimmed_gene.fastq \
    /app/data/gene_reference.fasta \
    /app/output/gene_analysis \
    8 \
    gene_name
```

### Batch Processing

```bash
docker run -it \
    -v /path/to/your/data:/app/data \
    -v /path/to/output:/app/output \
    cocoa-project \
    batch_devider \
    /app/data/gene_list.tsv \
    /app/data \
    /app/output/batch_results \
    /app/data/demultiplexed_reads \
    /app/data/reference_files \
    8
```
