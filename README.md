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
- **Multi-platform**: Supports both AMD64 and ARM64 architectures
- **Automated CI/CD**: Continuous integration and deployment via GitHub Actions
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
- For haplotype analysis: Docker socket access
- Sufficient disk space (recommend 50GB+ for full analysis)
- Internet connection for downloading the Docker image

### Pull the Docker Image

```bash  
# Pull the latest stable version  
docker pull ghcr.io/alienshe/cocoa_project_scripts:stable  

# Or pull a specific version  
docker pull ghcr.io/alienshe/cocoa_project_scripts:v1.0.0
