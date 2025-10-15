# POSEIDON Cancer Genomics Pipeline

A comprehensive cancer genomics data processing pipeline providing **two complementary approaches** for RNA-seq data analysis:

1. **Manual Pipeline** - Excel-based, manually curated sample selection
2. **Automated Pipeline** - Database-driven, automated discovery and classification

## Features

### Manual Pipeline
- **Excel-Based Curation**: Pre-selected samples organized by cancer type and category
- **High Control**: Manual review and selection of samples
- **Fast Processing**: Optimized for known, curated datasets
- **4 Sample Categories**: Tumors, Controls, Bulk_CellTypes, Premalignant

### Automated Pipeline
- **Maximum Metadata Collection**: Downloads ALL 192 available ENA fields
- **Multi-Source Integration**: Combines ENA, RunInfo, BioSample, BioProject, GEO, SRA XML, and ffq data
- **Cancer Classification**: Classifies samples into tumor/normal/cell line categories
- **Comprehensive Coverage**: 304 samples with 238+ columns of metadata

## Quick Start

### Manual Pipeline (Excel-Based)
```bash
# 1) Generate sample lists from Excel files
python scripts/manual_pipeline/GEO_sampleSetup_enhanced_VP.py data/manual_metadata/Pancreas.xlsx

# 2) Download SRA files
cd POSEIDON/Tumors/Pancreas/
./sratoolkit.sh | bsub

# 3) Convert to FASTQ
for i in *.sra; do ./fdump.sh $i | bsub; done

# 4) STAR alignment
for i in *1.fastq.gz; do bash star_2pass-paired.sh $i | bsub; done
```

### Automated Pipeline (Database-Driven)
```bash
# Test mode (safe limits)
python scripts/automated_pipeline/cancer_analysis_pipeline.py -c "lung cancer" --test

# Full run (be mindful of storage)
python scripts/automated_pipeline/cancer_analysis_pipeline.py -c "esophageal adenocarcinoma"
```

### Individual Automated Steps
```bash
# 1) Search SRR IDs by cancer type
python scripts/automated_pipeline/cancer_type_search.py -c "lung cancer" -o lung_srr_list.txt

# 2) Collect metadata (ENA + optional sources)
bash scripts/automated_pipeline/comprehensive_metadata_pipeline.sh -i lung_srr_list.txt -o output/lung_metadata --classify

# 3) Download FASTQs by category (optional)
python scripts/automated_pipeline/download_fastq_by_category.py -i output/lung_metadata/classified_metadata.tsv -o fastq_downloads
```

### Output
- **Manual Pipeline**: `POSEIDON/<category>/<cancer_type>/` directories with processed data
- **Automated Pipeline**: `output/<analysis>/comprehensive_metadata.tsv` – merged metadata
- **Both**: Ready for downstream analysis (junction analysis, AltAnalyze)

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

## Requirements

- Entrez Direct (esearch, efetch)
- curl, jq, python3
- pandas, numpy, requests (see `requirements.txt`)
- ffq (optional, for file information)

## Documentation

- **Comprehensive Guide**: `docs/COMPREHENSIVE_PIPELINE_GUIDE.md` - Complete documentation for both pipelines
- **Manual Pipeline**: Excel-based approach with pre-curated samples
- **Automated Pipeline**: Database-driven approach with comprehensive metadata collection

## Scripts

### Manual Pipeline
- `scripts/manual_pipeline/GEO_sampleSetup_enhanced_VP.py` – Generate sample lists from Excel
- `scripts/manual_pipeline/sratoolkit.sh` – Download SRA files
- `scripts/manual_pipeline/fdump.sh` – Convert SRA to FASTQ
- `scripts/manual_pipeline/star_2pass-paired.sh` – STAR 2-pass alignment

### Automated Pipeline
- `scripts/automated_pipeline/cancer_analysis_pipeline.py` – Orchestrates end-to-end workflow
- `scripts/automated_pipeline/cancer_type_search.py` – Finds SRR IDs by cancer type
- `scripts/automated_pipeline/comprehensive_metadata_pipeline.sh` – Collects and merges metadata
- `scripts/automated_pipeline/cancer_classification.py` – Classifies samples
- `scripts/automated_pipeline/download_fastq_by_category.py` – Downloads FASTQs by category
- `scripts/automated_pipeline/merge_metadata_maximum.py` – Merges all sources

## Data Sources

- **ENA**: 192 fields (maximum available)
- **RunInfo**: SRA run information
- **BioSample**: Clinical metadata (when available)
- **BioProject**: Study metadata (when available)
- **GEO**: Gene expression metadata (when available)
- **SRA XML**: Detailed technical metadata
- **ffq**: Comprehensive file information

## Documentation & Examples

- Guides: `docs/guides/` (installation, pipeline guide, setup)
- Reports: `docs/reports/` (test and validation reports)
- Examples and demo outputs: `examples/`
