# Complete Cancer Analysis Pipeline Guide

## Overview

This repository provides a comprehensive pipeline for cancer analysis that integrates four key steps:

1. **Cancer Type Search** - Find SRR IDs for specific cancer types
2. **Metadata Collection** - Gather comprehensive metadata from multiple sources
3. **Sample Classification** - Classify samples into categories (Tumor/Pre-malignant/Normal/Cell_line)
4. **FASTQ Download** - Download and organize FASTQ files by category

## Quick Start

### Complete Pipeline (Recommended)

```bash
# Run complete pipeline for a cancer type
python scripts/cancer_analysis_pipeline.py -c "lung cancer" --test

# Run with full results (be careful with storage!)
python scripts/cancer_analysis_pipeline.py -c "esophageal adenocarcinoma"
```

### Individual Steps

```bash
# Step 1: Search for cancer type
python scripts/cancer_type_search.py -c "lung cancer" -o lung_cancer_srr_list.txt

# Step 2: Collect metadata
bash scripts/comprehensive_metadata_pipeline.sh lung_cancer_srr_list.txt metadata.tsv

# Step 3: Classify samples
python scripts/cancer_classification.py -i metadata.tsv -o classified_metadata.tsv

# Step 4: Download FASTQ files
python scripts/download_fastq_by_category.py -i classified_metadata.tsv -o fastq_downloads
```

## Directory Structure

The pipeline creates an organized directory structure:

```
cancer_analysis_output/
└── {cancer_type}/
    ├── metadata/                    # SRR lists and comprehensive metadata
    │   ├── {cancer_type}_srr_list.txt
    │   └── comprehensive_metadata.tsv
    ├── classification/             # Classified metadata results
    │   └── classified_metadata.tsv
    ├── fastq_downloads/            # Organized FASTQ files by category
    │   ├── Tumor/
    │   │   ├── raw/                # Raw FASTQ files
    │   │   ├── processed/          # Processed files
    │   │   ├── logs/               # Download logs
    │   │   └── tumor_srr_list.txt  # SRR list for this category
    │   ├── Pre-malignant/
    │   ├── Normal/
    │   ├── Cell_line/
    │   └── Unknown/
    └── logs/                       # Pipeline execution logs
```

## Scripts Overview

### 1. Cancer Type Search (`cancer_type_search.py`)

Searches SRA/ENA databases for cancer types and extracts associated SRR IDs.

**Usage:**
```bash
python scripts/cancer_type_search.py -c "lung cancer" -o output.txt -m 1000
```

**Parameters:**
- `-c, --cancer-type`: Cancer type to search for
- `-o, --output`: Output file for SRR IDs
- `-m, --max-results`: Maximum results per source (default: 1000)
- `--test`: Run in test mode with limited results

### 2. Metadata Collection (`comprehensive_metadata_pipeline.sh`)

Collects comprehensive metadata from multiple sources (ENA, NCBI, GEO, SRA).

**Usage:**
```bash
bash scripts/comprehensive_metadata_pipeline.sh srr_list.txt output_metadata.tsv
```

### 3. Cancer Classification (`cancer_classification.py`)

Classifies samples into categories using comprehensive metadata and cell line reference database.

**Usage:**
```bash
python scripts/cancer_classification.py -i metadata.tsv -o classified_metadata.tsv
```

**Classification Categories:**
- **Tumor**: Cancer/malignant samples
- **Pre-malignant**: Barrett's esophagus, dysplasia, etc.
- **Normal**: Healthy/control samples
- **Cell_line**: Cell culture samples
- **Unknown**: Unclassified samples

### 4. FASTQ Download (`download_fastq_by_category.py`)

Downloads FASTQ files organized by classification categories.

**Usage:**
```bash
python scripts/download_fastq_by_category.py -i classified_metadata.tsv -o fastq_downloads -w 4
```

**Parameters:**
- `-i, --input`: Input classified metadata TSV file
- `-o, --output-dir`: Output directory for FASTQ files
- `-w, --workers`: Number of parallel download workers (default: 4)
- `--max-per-category`: Maximum downloads per category (for testing)
- `--categories`: Categories to download (default: all except Unknown)

### 5. Complete Pipeline (`cancer_analysis_pipeline.py`)

Orchestrates all four steps in a single command.

**Usage:**
```bash
python scripts/cancer_analysis_pipeline.py -c "lung cancer" -o output_dir --test
```

**Parameters:**
- `-c, --cancer-type`: Cancer type to analyze
- `-o, --output-dir`: Output directory (default: cancer_analysis_output)
- `-m, --max-results`: Maximum SRR results to search (default: 1000)
- `--max-per-category`: Maximum FASTQ downloads per category
- `--skip-download`: Skip FASTQ download step
- `--test`: Run in test mode with limited results

## Installation & Setup

### Prerequisites

1. **Python 3.8+** with virtual environment
2. **SRA Toolkit** for FASTQ downloads
3. **Required Python packages**

### Setup

```bash
# Clone repository
git clone <repository-url>
cd comprehensive-sra-metadata-pipeline

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install Python dependencies
pip install -r requirements.txt

# Install SRA Toolkit (macOS)
brew install sratoolkit

# Install SRA Toolkit (Linux)
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
export PATH=$PATH:$(pwd)/sratoolkit.current-ubuntu64/bin
```

## Examples

### Example 1: Test Run with Lung Cancer

```bash
# Run complete pipeline in test mode
python scripts/cancer_analysis_pipeline.py -c "lung cancer" --test

# This will:
# 1. Search for lung cancer samples (limited to 50)
# 2. Collect metadata for found samples
# 3. Classify samples into categories
# 4. Download FASTQ files (max 5 per category)
```

### Example 2: Esophageal Cancer Analysis

```bash
# Run complete pipeline for esophageal cancer
python scripts/cancer_analysis_pipeline.py -c "esophageal adenocarcinoma" -o esophageal_analysis

# Check results
ls esophageal_analysis/esophageal_adenocarcinoma/
```

### Example 3: Skip FASTQ Download

```bash
# Run pipeline without downloading FASTQ files
python scripts/cancer_analysis_pipeline.py -c "breast cancer" --skip-download
```

## Important Notes

### dbGaP Authorization

Some SRR IDs require dbGaP authorization. The pipeline handles this gracefully:
- Samples requiring authorization are skipped with a warning
- Download continues with available public samples
- Authorization status is logged for each sample

### Storage Requirements

FASTQ files are large. Consider:
- **Test mode**: Use `--test` flag for limited downloads
- **Category limits**: Use `--max-per-category` to limit downloads
- **Storage space**: Ensure sufficient disk space (GB to TB depending on sample count)

### Performance

- **Parallel downloads**: Adjust `-w` parameter based on system capabilities
- **Network limits**: SRA has rate limits; respect them
- **Resume capability**: Downloads can be resumed if interrupted

## Troubleshooting

### Common Issues

1. **SRA Toolkit not found**
   ```bash
   # Install SRA Toolkit
   brew install sratoolkit  # macOS
   # or download from NCBI website
   ```

2. **Permission denied errors**
   - Some samples require dbGaP authorization
   - Check logs for authorization requirements
   - Apply for dbGaP access if needed

3. **Out of disk space**
   - Use `--max-per-category` to limit downloads
   - Monitor disk usage during downloads
   - Consider downloading to external storage

4. **Network timeouts**
   - Reduce parallel workers (`-w` parameter)
   - Check network connection
   - Retry failed downloads

### Logs

Check log files for detailed information:
- `fastq_download.log`: Download progress and errors
- `cancer_analysis_pipeline.log`: Pipeline execution logs
- Category-specific logs in `fastq_downloads/{category}/logs/`

## Output Files

### Metadata Files
- `comprehensive_metadata.tsv`: Raw metadata from all sources
- `classified_metadata.tsv`: Metadata with classification results

### SRR Lists
- `{category}_srr_list.txt`: SRR IDs for each category
- Useful for reproducibility and downstream analysis

### Summary Reports
- `pipeline_summary.txt`: Complete pipeline results
- `download_summary.txt`: FASTQ download statistics

## Contributing

When adding new features:
1. Test with small datasets first
2. Add comprehensive error handling
3. Update documentation
4. Ensure compatibility with existing pipeline

## License

[Add your license information here]