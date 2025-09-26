# Cancer Type-Based Metadata Collection System

This branch contains an enhanced version of the metadata collection system that allows users to search for cancer types by name instead of requiring pre-existing SRR IDs.

## Overview

The new system enables researchers to:
1. **Search by cancer type name** (e.g., "esophageal adenocarcinoma", "lung squamous cell carcinoma")
2. **Automatically discover relevant studies** in SRA/ENA databases
3. **Extract SRR IDs** from matching studies
4. **Collect comprehensive metadata** from multiple sources
5. **Apply cancer classification** (tumor/normal/cell line)

## Key Features

### 🎯 Cancer Type Search
- Accepts natural language cancer type descriptions
- Supports multiple cancer types and synonyms
- Searches across SRA and ENA databases
- Handles variations in terminology (cancer vs carcinoma, etc.)

### 📊 Comprehensive Metadata Collection
- **SRA RunInfo**: Basic run information
- **ENA Filereport**: Enhanced metadata with cancer-specific fields
- **BioSample**: Clinical and sample metadata
- **BioProject**: Study-level information
- **GEO**: Gene Expression Omnibus metadata (optional)
- **SRA XML**: Detailed XML metadata (optional)
- **ffq**: Fastq file information (optional)

### 🏷️ Cancer Classification
- **Tumor**: Malignant samples
- **Normal**: Healthy/control samples
- **Cell line**: Cultured cell samples
- **Pre-malignant**: Precancerous conditions (e.g., Barrett's esophagus)
- **Unknown**: Unclassified samples

## Usage

### Basic Usage

```bash
# Search for esophageal adenocarcinoma
scripts/collect_cancer_metadata.sh -c "esophageal adenocarcinoma"

# Search for lung squamous cell carcinoma with all metadata sources
scripts/collect_cancer_metadata.sh -c "lung squamous cell carcinoma" --with-geo --with-xml --with-ffq

# Test mode with limited results
scripts/collect_cancer_metadata.sh -c "breast cancer" --test
```

### Advanced Options

```bash
scripts/collect_cancer_metadata.sh -c "pancreatic cancer" \
  --output-dir ./output/pancreatic_study \
  --max-results 2000 \
  --with-geo \
  --with-xml \
  --with-ffq
```

### Command Line Options

| Option | Description |
|--------|-------------|
| `-c, --cancer-type` | **Required**: Cancer type to search for |
| `-o, --output-dir` | Output directory (default: ./output/cancer_type_collection) |
| `--with-geo` | Include GEO metadata collection |
| `--with-xml` | Include SRA XML metadata collection |
| `--with-ffq` | Include ffq metadata collection |
| `--max-results` | Maximum results per source (default: 1000) |
| `--test` | Run in test mode (limited results) |
| `--skip-search` | Skip cancer type search (use existing SRR list) |
| `--skip-merge` | Skip metadata merging step |
| `--skip-classify` | Skip cancer classification step |

## Workflow

The system follows a 4-step pipeline:

### Step 1: Cancer Type Search
- Searches SRA/ENA databases for studies matching the cancer type
- Extracts SRR IDs from matching studies
- Generates `cancer_type_srr_list.txt`

### Step 2: Metadata Collection
- Collects metadata from multiple sources for each SRR ID
- Saves raw data in `raw/` directory
- Sources: RunInfo, ENA, BioSample, BioProject, GEO, SRA XML, ffq

### Step 3: Metadata Merging
- Merges all metadata sources into unified dataset
- Creates `ultimate_metadata.tsv`
- Handles duplicate columns and missing data

### Step 4: Cancer Classification
- Applies cancer classification algorithm
- Creates `classified_metadata.tsv`
- Provides summary statistics

## Output Files

```
output/cancer_type_collection/
├── cancer_type_srr_list.txt     # List of SRR IDs found
├── ultimate_metadata.tsv         # Merged metadata from all sources
├── classified_metadata.tsv      # Metadata with cancer classifications
└── raw/                         # Raw metadata from each source
    ├── runinfo.csv
    ├── ena_read_run.tsv
    ├── biosample.jsonl
    ├── bioproject.jsonl
    ├── geo_metadata.tsv
    ├── sra_xml.jsonl
    └── ffq.jsonl
```

## Supported Cancer Types

The system supports various cancer types including:

- **Esophageal**: esophageal adenocarcinoma, esophageal squamous cell carcinoma
- **Lung**: lung adenocarcinoma, lung squamous cell carcinoma, small cell lung cancer
- **Breast**: breast cancer, breast carcinoma, triple-negative breast cancer
- **Colorectal**: colorectal cancer, colon cancer, rectal cancer
- **Prostate**: prostate cancer, prostate adenocarcinoma
- **Pancreatic**: pancreatic cancer, pancreatic adenocarcinoma
- **Liver**: hepatocellular carcinoma, liver cancer
- **Kidney**: renal cell carcinoma, kidney cancer
- **Brain**: glioma, glioblastoma, brain cancer
- **Skin**: melanoma, cutaneous melanoma
- **Blood**: leukemia, lymphoma

## Cancer Classification Fields

The system adds the following classification fields:

| Field | Description | Values |
|-------|-------------|---------|
| `top_label` | Primary classification | Tumor, Normal, Cell line, Pre-malignant, Unknown |
| `is_cell_line` | Is this a cell line? | yes, no |
| `is_bulk_sorted` | Is this bulk/sorted? | yes, no |
| `is_control` | Is this a control sample? | yes, no |
| `adjacent_normal` | Is this adjacent normal tissue? | yes, no |
| `barrett_grade` | Barrett's esophagus grade | LGD, HGD, indefinite, no dysplasia, unknown |
| `tissue_origin` | Tissue of origin | esophagus, stomach, lung, breast, etc. |

## Examples

### Example 1: Esophageal Adenocarcinoma
```bash
scripts/collect_cancer_metadata.sh -c "esophageal adenocarcinoma" --test
```

**Results:**
- Found 50 SRR IDs
- All classified as "Tumor"
- Tissue origin: "esophagus"
- Study: "Shifts in Serum Bile Acid Profiles Associated with Barretts Esophagus and Stages of Progression to Esophageal Adenocarcinoma"

### Example 2: Lung Cancer with Full Metadata
```bash
scripts/collect_cancer_metadata.sh -c "lung squamous cell carcinoma" \
  --with-geo --with-xml --with-ffq \
  --max-results 500
```

**Results:**
- Comprehensive metadata collection
- GEO annotations included
- SRA XML details included
- Fastq file information included

## Comparison with Original System

| Feature | Original System | New System |
|---------|----------------|------------|
| **Input** | SRR IDs (pre-existing list) | Cancer type names |
| **Discovery** | Manual curation required | Automatic database search |
| **Flexibility** | Limited to known SRR IDs | Discovers new studies |
| **Usability** | Requires bioinformatics knowledge | User-friendly cancer terminology |
| **Coverage** | Limited to curated datasets | Broader database coverage |

## Technical Details

### Search Strategy
1. **Term Normalization**: Converts cancer types to multiple search variations
2. **Database Queries**: Searches SRA and ENA using multiple strategies
3. **Result Aggregation**: Combines results from multiple sources
4. **Deduplication**: Removes duplicate SRR IDs

### Metadata Sources
- **SRA**: Sequence Read Archive (NCBI)
- **ENA**: European Nucleotide Archive (EBI)
- **GEO**: Gene Expression Omnibus (NCBI)
- **BioSample**: Sample metadata (NCBI)
- **BioProject**: Project metadata (NCBI)

### Error Handling
- Graceful handling of API failures
- Rate limiting to respect database limits
- Fallback strategies for different data sources
- Comprehensive logging and error reporting

## Requirements

- Python 3.6+
- Bash shell
- Internet connection for database queries
- Optional: EDirect tools for enhanced SRA access

## Installation

1. Clone the repository
2. Ensure Python 3.6+ is available
3. Make scripts executable: `chmod +x scripts/*.sh`
4. Test with: `scripts/collect_cancer_metadata.sh -c "esophageal adenocarcinoma" --test`

## Troubleshooting

### Common Issues

1. **No results found**: Try broader search terms or use `--broad` flag
2. **API errors**: Check internet connection and try again
3. **Permission errors**: Ensure scripts are executable
4. **Python errors**: Verify Python 3.6+ is available

### Debug Mode

Use `--test` flag for limited results and easier debugging.

## Future Enhancements

- [ ] Support for additional cancer types
- [ ] Integration with TCGA and other cancer databases
- [ ] Machine learning-based cancer classification
- [ ] Web interface for easier access
- [ ] Batch processing for multiple cancer types
- [ ] Integration with clinical trial databases

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this system in your research, please cite:

```
Cancer Type-Based Metadata Collection System
Valerii Try, 2024
https://github.com/your-repo/cancer-metadata-collection
```

## Support

For questions or issues, please:
1. Check the troubleshooting section
2. Review the examples
3. Open an issue on GitHub
4. Contact the maintainers