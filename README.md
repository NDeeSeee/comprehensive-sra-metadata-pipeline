# Comprehensive SRA Metadata Pipeline

A comprehensive pipeline for downloading and processing SRA metadata with maximum field coverage and cancer classification capabilities.

## Features

- **Maximum Metadata Collection**: Downloads ALL 192 available ENA fields
- **Multi-Source Integration**: Combines ENA, RunInfo, BioSample, BioProject, GEO, SRA XML, and ffq data
- **Cancer Classification**: Classifies samples into tumor/normal/cell line categories
- **Comprehensive Coverage**: 304 samples with 238+ columns of metadata

## Quick Start

### Basic Usage
```bash
# Download comprehensive metadata
./scripts/comprehensive_metadata_pipeline.sh -i data/srr_list.txt

# Download metadata + cancer classification
./scripts/comprehensive_metadata_pipeline.sh -i data/srr_list.txt --classify
```

### Output Files
- `results/comprehensive_metadata.tsv` - Maximum metadata (238 columns)
- `results/classified_metadata.tsv` - Metadata + cancer classification (99 columns)

## Requirements

- Entrez Direct (esearch, efetch)
- curl, jq, python3
- pandas (included in venv/)
- ffq (optional, for file information)

## Scripts

- `scripts/comprehensive_metadata_pipeline.sh` - Main pipeline script
- `scripts/cancer_classification_enhanced.py` - Cancer classification script
- `scripts/merge_metadata_maximum.py` - Metadata merging script

## Data Sources

- **ENA**: 192 fields (maximum available)
- **RunInfo**: SRA run information
- **BioSample**: Clinical metadata (when available)
- **BioProject**: Study metadata (when available)
- **GEO**: Gene expression metadata (when available)
- **SRA XML**: Detailed technical metadata
- **ffq**: Comprehensive file information

## Results

Current dataset: **304 samples** with **238 columns** of comprehensive metadata.

Cancer classification results:
- 293 tumors
- 11 cell lines
- 24 adjacent normal samples
- All samples correctly identified as esophagus origin
