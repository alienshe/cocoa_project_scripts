# Cocoa Project Genomics Pipeline

A comprehensive Docker-based bioinformatics pipeline for demultiplexing nanopore sequencing data and performing haplotype analysis on plant genomic data. This pipeline was developed for analyzing NLR (Nucleotide-binding Leucine-rich Repeat) genes in cocoa (*Theobroma cacao*) but can be adapted for other genomic studies.

## Overview

The Cocoa Project Genomics Pipeline provides a complete workflow for:

1. **Demultiplexing**: Separating pooled nanopore sequencing reads based on barcodes
2. **Quality Control**: Adapter trimming and read filtering
3. **Alignment**: Mapping reads to reference sequences using minimap2
4. **Variant Calling**: Identifying genetic variants using Clair3
5. **Haplotype Reconstruction**: Reconstructing haplotypes using devider
6. **Batch Processing**: Analyzing multiple samples simultaneously

## Key Features

- **Containerized**: Complete pipeline packaged in Docker for reproducibility
- **Real Demo Data**: Includes actual pangkep study data for testing and learning
- **Multi-platform**: Supports both AMD64 and ARM64 architectures (with platform flag for ARM64)
- **Self-contained**: All tools and models included in the container
- **Interactive Menu**: User-friendly command interface with built-in help
- **Comprehensive Output**: Detailed logs, alignments, variants, and haplotypes
- **Batch Processing**: Efficient analysis of multiple samples

## Scientific Background

This pipeline was designed for studying NLR genes, which are important components of plant immune systems. The workflow processes long-read nanopore sequencing data to:

- Identify genetic variants in NLR genes
- Reconstruct haplotypes to understand allelic diversity
- Analyze structural variations that may affect disease resistance
- Support breeding programs and evolutionary studies

## Installation

### Prerequisites

- Docker installed on your system
- Sufficient disk space (recommend 50GB+ for full analysis)
- Internet connection for downloading the Docker image

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

# For x86_64 systems (Linux/Intel Mac):
docker build -t cocoa-pipeline .

# For ARM64 systems (Apple Silicon Mac):
docker build --platform=linux/amd64 -t cocoa-pipeline .
```

## Quick Start

### Interactive Menu

```bash
# Run the container to see the interactive menu
docker run -it cocoa-pipeline

# Or with volume mounting for persistent data
docker run -it -v $(pwd)/data:/app/data cocoa-pipeline
```

The interactive menu provides:

Demo commands: Try the pipeline with included test data
- Production commands: Run on your own data
- Tool access: Direct access to all bioinformatics tools
- Status checks: Verify all tools and models are properly installed

### Demo Analysis

```bash
# List available demo files
docker run -it cocoa-pipeline list_demo_files

# Run demultiplexing demo
docker run -it -v $(pwd)/output:/app/output cocoa-pipeline \
    demo_demultiplex /app/output/demux_demo

# Run single gene analysis demo
docker run -it -v $(pwd)/output:/app/output cocoa-pipeline \
    demo_devider Tc1318 /app/output/gene_demo

# Run batch processing demo
docker run -it -v $(pwd)/output:/app/output cocoa-pipeline \
    demo_batch_devider /app/output/batch_demo
```

## Usage

### Demultiplexing

```bash

docker run -it \
    -v /path/to/your/data:/app/data \
    -v /path/to/output:/app/output \
    cocoa-pipeline \
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
    cocoa-pipeline \
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
    cocoa-pipeline \
    batch_devider \
    /app/data/gene_list.tsv \
    /app/data \
    /app/output/batch_results \
    /app/data/demultiplexed_reads \
    /app/data/reference_files \
    8
```

### Interactive Shell

```bash
# Access the container for custom analysis
docker run -it -v $(pwd)/data:/app/data cocoa-pipeline bash

# Inside the container, all tools are available:
# porechop, minimap2, samtools, devider, clair3
```

## Tools Included

The pipeline includes all necessary bioinformatics tools:

porechop: Adapter trimming for nanopore reads
minimap2: Fast alignment of long reads
samtools: SAM/BAM file manipulation
tabix: Indexing for genomic data
lofreq: Variant calling
devider: Haplotype reconstruction
clair3: Deep learning-based variant calling (with r1041_e82_400bps_hac_v500 model)

## Architecture

The pipeline is built on:

Base: Ubuntu Linux via micromamba container
Package Manager: Micromamba for fast, reliable package management
Sources: Bioconda and conda-forge for bioinformatics tools
Models: Pre-installed Clair3 model for nanopore variant calling

## File Structure

```bash
cocoa_project_scripts/
├── demo_files/                    # Test data and examples
│   ├── inputs/
│   │   ├── pangkep_raw_reads/
│   │   ├── pangkep_run_bc01_demultiplexed/
│   │   └── reference_files/
│   └── outputs/
├── scripts/
│   ├── full_demultiplex_inc_minimap.sh
│   ├── sylvies_devider_pipeline.sh
│   └── batch_run_devider_pipeline.sh
├── Dockerfile
└── README.md
```

### Platform Support

```bash
Linux (x86_64): Native support
macOS (Intel): Native support
macOS (Apple Silicon): Use --platform=linux/amd64 flag
Windows: Via Docker Desktop
```

### Performance Notes

```bash
Memory: Minimum 8GB RAM recommended, 16GB+ for large datasets
CPU: Multi-threading supported (adjust thread count in commands)
Storage: SSD recommended for faster I/O operations
Network: Initial download ~2-3GB for container image
```

### Troubleshooting

Common Issues
Permission errors: Ensure Docker has access to mounted directories
Memory issues: Increase Docker memory limits for large datasets
Platform issues: Use --platform=linux/amd64 on Apple Silicon
Tool not found: Check tool availability with the interactive menu

### Getting Help

```bash
# Check tool status
docker run -it cocoa-pipeline

# Access interactive shell for debugging
docker run -it cocoa-pipeline bash

# Check individual tool help
docker run -it cocoa-pipeline bash -c "devider --help"
```

