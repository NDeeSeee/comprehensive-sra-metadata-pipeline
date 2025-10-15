# POSEIDON Comprehensive Pipeline Guide

## Overview

The POSEIDON project provides **two complementary approaches** for cancer genomics data processing:

1. **Manual Pipeline** - Excel-based, manually curated sample selection
2. **Automated Pipeline** - Database-driven, automated discovery and classification

Both pipelines ultimately produce the same output: processed RNA-seq data ready for downstream analysis.

## Directory Structure

```
ValeriiGitRepo/
├── scripts/
│   ├── automated_pipeline/          # Database-driven approach
│   │   ├── cancer_analysis_pipeline.py
│   │   ├── cancer_classification.py
│   │   ├── cancer_type_search.py
│   │   ├── comprehensive_metadata_pipeline.sh
│   │   ├── download_fastq_by_category.py
│   │   ├── extract_zero_tcga_cancers.py
│   │   └── merge_metadata_maximum.py
│   ├── manual_pipeline/             # Excel-based approach
│   │   ├── GEO_sampleSetup_enhanced_VP.py
│   │   ├── sratoolkit.sh
│   │   ├── fdump.sh
│   │   └── star_2pass-paired.sh
│   └── downstream_analysis/         # Post-alignment analysis
├── data/
│   ├── automated_metadata/         # Auto-generated metadata
│   └── manual_metadata/            # Manually curated Excel files
│       ├── Esophagus.xlsx
│       ├── Oropharyngeal.xlsx
│       ├── Gallbladder.xlsx
│       └── ... (all cancer types)
└── docs/
    ├── COMPREHENSIVE_PIPELINE_GUIDE.md
    └── README.md
```

## Pipeline Comparison

| Feature | Manual Pipeline | Automated Pipeline |
|---------|----------------|-------------------|
| **Sample Selection** | Manual curation in Excel | Automated database search |
| **Metadata Source** | Pre-selected studies | Comprehensive ENA/SRA search |
| **Control** | High (manual review) | Medium (automated classification) |
| **Speed** | Fast (pre-selected) | Slower (comprehensive search) |
| **Coverage** | Limited to known studies | Maximum possible coverage |
| **Maintenance** | High (manual updates) | Low (automated) |
| **Best For** | Focused studies, known samples | Discovery, comprehensive analysis |

## Manual Pipeline (Excel-Based)

### Overview
The manual pipeline uses **pre-curated Excel files** containing manually selected samples organized by cancer type and sample category.

### Workflow

#### 1. **Metadata Preparation**
- **Input**: Excel files in `data/manual_metadata/`
- **Structure**: 4 sheets per cancer type:
  - `Tumors`: Malignant samples
  - `Controls`: Healthy/control samples
  - `Bulk_CellTypes`: Bulk cell type samples
  - `Premalignant`: Pre-malignant samples

#### 2. **Sample List Generation**
```bash
python scripts/manual_pipeline/GEO_sampleSetup_enhanced_VP.py data/manual_metadata/Pancreas.xlsx
```
- **Output**: `sample_list.txt` files in POSEIDON directories
- **Format**: `BioSample_ID\tFASTQ_R1\tFASTQ_R2`

#### 3. **SRA Download**
```bash
cd POSEIDON/Tumors/Pancreas/
./sratoolkit.sh | bsub
```
- **Tool**: SRA Toolkit `prefetch`
- **Input**: `sample_list.txt` or SRA list file
- **Output**: `.sra` files

#### 4. **FASTQ Conversion**
```bash
for i in *.sra; do ./fdump.sh $i | bsub; done
```
- **Tool**: SRA Toolkit `fastq-dump`
- **Input**: `.sra` files
- **Output**: `*_1.fastq.gz` and `*_2.fastq.gz` files

#### 5. **STAR Alignment**
```bash
for i in *1.fastq.gz; do bash star_2pass-paired.sh $i | bsub; done
```
- **Tool**: STAR 2-pass alignment
- **Reference**: GRCh38
- **Input**: FASTQ files
- **Output**: BAM files

### Manual Pipeline Scripts

#### `GEO_sampleSetup_enhanced_VP.py`
- **Purpose**: Generate sample_list.txt from Excel metadata
- **Input**: Excel file with 4 sheets
- **Output**: sample_list.txt files in POSEIDON directories
- **Usage**: `python GEO_sampleSetup_enhanced_VP.py metadata_file.xlsx`

#### `sratoolkit.sh`
- **Purpose**: Download SRA files using prefetch
- **Input**: SRA list file
- **Output**: .sra files
- **Usage**: `./sratoolkit.sh | bsub`

#### `fdump.sh`
- **Purpose**: Convert SRA to FASTQ
- **Input**: .sra files
- **Output**: *_1.fastq.gz, *_2.fastq.gz files
- **Usage**: `for i in *.sra; do ./fdump.sh $i | bsub; done`

#### `star_2pass-paired.sh`
- **Purpose**: STAR 2-pass alignment
- **Input**: FASTQ files
- **Output**: BAM files
- **Usage**: `for i in *1.fastq.gz; do bash star_2pass-paired.sh $i | bsub; done`

## Automated Pipeline (Database-Driven)

### Overview
The automated pipeline uses **comprehensive database searches** to discover and classify samples automatically.

### Workflow

#### 1. **Cancer Type Search**
```bash
python scripts/automated_pipeline/cancer_type_search.py -c "pancreatic cancer" -o output.txt
```
- **Purpose**: Search SRA/ENA databases by cancer type
- **Input**: Cancer type name
- **Output**: List of SRR IDs

#### 2. **Comprehensive Metadata Collection**
```bash
bash scripts/automated_pipeline/comprehensive_metadata_pipeline.sh
```
- **Sources**: ENA (192 fields), RunInfo, BioSample, BioProject, GEO, SRA XML, ffq
- **Output**: `comprehensive_metadata.tsv` (238+ columns)

#### 3. **Metadata Merging**
```bash
python scripts/automated_pipeline/merge_metadata_maximum.py -i raw/ -o ultimate_metadata.tsv
```
- **Purpose**: Merge all metadata sources
- **Input**: Raw metadata files
- **Output**: Comprehensive metadata TSV

#### 4. **Cancer Classification**
```bash
python scripts/automated_pipeline/cancer_classification.py -i metadata.tsv -o classified.tsv
```
- **Purpose**: Classify samples into tumor/normal/cell line categories
- **Features**: Cell line reference database integration
- **Output**: `classified_metadata.tsv`

#### 5. **FASTQ Download by Category**
```bash
python scripts/automated_pipeline/download_fastq_by_category.py -i classified.tsv -o fastq_downloads/
```
- **Purpose**: Download FASTQ files organized by classification
- **Input**: Classified metadata
- **Output**: Organized FASTQ files

### Automated Pipeline Scripts

#### `cancer_type_search.py`
- **Purpose**: Search SRA/ENA for cancer types
- **Input**: Cancer type name
- **Output**: SRR ID list
- **Usage**: `python cancer_type_search.py -c "pancreatic cancer" -o output.txt`

#### `comprehensive_metadata_pipeline.sh`
- **Purpose**: Collect metadata from all sources
- **Sources**: ENA, RunInfo, BioSample, BioProject, GEO, SRA XML, ffq
- **Output**: Comprehensive metadata files

#### `merge_metadata_maximum.py`
- **Purpose**: Merge all metadata sources
- **Input**: Raw metadata files
- **Output**: Comprehensive metadata TSV
- **Usage**: `python merge_metadata_maximum.py -i raw/ -o ultimate_metadata.tsv`

#### `cancer_classification.py`
- **Purpose**: Classify samples by type
- **Input**: Comprehensive metadata
- **Output**: Classified metadata
- **Usage**: `python cancer_classification.py -i metadata.tsv -o classified.tsv`

#### `download_fastq_by_category.py`
- **Purpose**: Download FASTQ by classification
- **Input**: Classified metadata
- **Output**: Organized FASTQ files
- **Usage**: `python download_fastq_by_category.py -i classified.tsv -o fastq_downloads/`

## Downstream Analysis

### Junction Analysis
```bash
/software/LabShellScripts/RunAltAnalyze.from.bams/index-junction_hg38.sh
```

### AltAnalyze Pipeline
```bash
/software/LabShellScripts/RunAltAnalyze.from.bams/AltAnalyze.sh
```

## Technical Specifications

### Software Versions
- **SRA Toolkit**: 3.1.1 (prefetch), 2.10.4 (fastq-dump)
- **STAR**: 2.4.0h (2-pass alignment)
- **Reference Genome**: GRCh38
- **Compute Resources**: LSF job scheduler, 4-8 cores, 10-128GB RAM

### Data Storage
- **Organization**: By cancer type and sample category
- **File Formats**: SRA, FASTQ, BAM, Excel metadata
- **Storage**: Organized by cancer type and sample category

## Usage Recommendations

### Choose Manual Pipeline When:
- You have specific studies/samples in mind
- You need high control over sample selection
- You want to process known, curated datasets
- You need fast processing of pre-selected samples

### Choose Automated Pipeline When:
- You want maximum coverage of available data
- You're doing discovery-based research
- You want to minimize manual curation effort
- You need comprehensive metadata collection

### Hybrid Approach:
- Use automated pipeline for discovery
- Use manual pipeline for focused analysis
- Combine results from both approaches

## Getting Started

### Manual Pipeline Quick Start:
1. Prepare Excel files in `data/manual_metadata/`
2. Run `GEO_sampleSetup_enhanced_VP.py` for each cancer type
3. Execute download and alignment scripts
4. Proceed to downstream analysis

### Automated Pipeline Quick Start:
1. Run `cancer_type_search.py` for your cancer type
2. Execute `comprehensive_metadata_pipeline.sh`
3. Run `merge_metadata_maximum.py` and `cancer_classification.py`
4. Use `download_fastq_by_category.py` for organized downloads
5. Proceed to downstream analysis

## Support and Documentation

- **Main README**: `README.md`
- **Changelog**: `data/manual_metadata/CHANGELOG.md`
- **Examples**: `examples/` directory
- **Issues**: Check project documentation for troubleshooting

---

*This comprehensive guide covers both manual and automated approaches to cancer genomics data processing in the POSEIDON project.*