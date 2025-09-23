# Comprehensive SRA Metadata Collection Pipeline

A robust, production-ready pipeline for collecting comprehensive metadata from multiple authoritative sources for SRA (Sequence Read Archive) samples.

## 🎯 Overview

This pipeline collects **all possible metadata** for SRA samples by integrating data from:
- **NCBI SRA** (via EDirect tools)
- **ENA** (European Nucleotide Archive)
- **BioSample** database
- **BioProject** database
- **GEO** (Gene Expression Omnibus) - when applicable
- **ffq** tool for additional JSON metadata

## 📁 Project Structure

```
├── scripts/
│   ├── build_metadata_enhanced.sh    # Main metadata collection script
│   └── merge_metadata_enhanced.py    # Metadata merging and processing script
├── data/
│   └── srr_list.txt                  # Input: List of SRR accessions (one per line)
├── output/
│   └── neo_meta_final/               # Final comprehensive metadata output
│       ├── raw/                       # Raw data from all sources
│       ├── logs/                      # Execution logs
│       └── ultimate_metadata_corrected.tsv  # Final merged metadata table
└── docs/                             # Documentation
```

## 🚀 Quick Start

### Prerequisites

1. **NCBI EDirect Tools**: Install and configure
   ```bash
   # Install EDirect
   sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
   
   # Configure API key (optional but recommended)
   export NCBI_API_KEY="your_api_key_here"
   ```

2. **Python Environment**: Python 3.7+ with required packages
   ```bash
   pip install pandas numpy requests
   ```

3. **System Tools**: `curl`, `jq`, `awk`, `sort`, `uniq`

### Usage

1. **Prepare input file**: Create `data/srr_list.txt` with one SRR accession per line
   ```
   SRR1234567
   SRR1234568
   SRR1234569
   ```

2. **Run comprehensive metadata collection**:
   ```bash
   cd scripts
   ./build_metadata_enhanced.sh -i ../data/srr_list.txt -o ../output/comprehensive_metadata --with-geo --with-xml --with-ffq
   ```

3. **Process and merge metadata**:
   ```bash
   python3 merge_metadata_enhanced.py ../output/comprehensive_metadata
   ```

## 📊 Output

The final output (`ultimate_metadata_corrected.tsv`) contains **92 comprehensive metadata fields** including:

### Core Identifiers
- `run_accession` - SRA run ID
- `experiment_accession` - SRA experiment ID  
- `sample_accession` - SRA sample ID
- `study_accession` - SRA study ID
- `BioSample` - BioSample database ID
- `BioProject` - BioProject database ID

### Technical Details
- `LibraryLayout` - Single/Paired end
- `LibraryStrategy` - RNA-Seq, WGS, etc.
- `Platform` - Sequencing platform
- `Model` - Instrument model
- `read_count`, `base_count` - Sequencing statistics

### Sample Characteristics
- `age`, `sex`, `disease`, `tissue_type`
- `cell_line`, `cell_type`, `strain`
- `host`, `host_body_site`, `host_genotype`
- `environment_biome`, `location`, `temperature`

### File Access
- `submitted_ftp`, `fastq_ftp` - Download URLs
- `download_path` - Direct file paths

### Publication Data
- `Study_Pubmed_id` - PubMed references
- `study_title`, `sample_title` - Descriptive titles
- `first_public`, `last_updated` - Publication dates

## 🔧 Advanced Features

### Command Line Options

```bash
./build_metadata_enhanced.sh [OPTIONS]

Options:
  -i FILE              Input file with SRR accessions (required)
  -o DIR               Output directory (default: ./meta_out)
  --with-geo           Include GEO metadata collection
  --with-xml           Include SRA XML format data
  --with-ffq           Include ffq JSON metadata
```

### Environment Configuration

The script automatically handles:
- **Proxy settings**: Automatically unsets proxy variables
- **Conda environments**: Uses system Python for reliability
- **Error handling**: Robust error recovery and logging
- **Rate limiting**: Built-in delays for API compliance

## 🛠️ Technical Details

### Data Sources Integration

1. **SRA RunInfo**: Basic run information via EDirect
2. **ENA Filereport**: Comprehensive metadata via ENA API
3. **BioSample**: Sample characteristics and attributes
4. **BioProject**: Project-level metadata
5. **SRA XML**: Detailed technical specifications
6. **GEO**: Gene expression study metadata (when applicable)
7. **ffq**: Additional JSON-formatted metadata

### Error Handling

- **Robust API calls**: Automatic retry with exponential backoff
- **Graceful degradation**: Continues processing if individual sources fail
- **Comprehensive logging**: Detailed logs for debugging
- **Data validation**: Ensures data integrity throughout pipeline

### Performance Optimizations

- **Parallel processing**: Concurrent API calls where possible
- **Caching**: Avoids redundant API requests
- **Memory efficient**: Streams large datasets
- **Timeout handling**: Prevents hanging processes

## 📈 Quality Assurance

### Data Validation
- **Completeness checks**: Verifies all samples processed
- **Cross-reference validation**: Ensures data consistency across sources
- **Duplicate detection**: Handles and reports duplicate entries
- **Format validation**: Ensures proper data types and formats

### Testing
- **Sample validation**: Tested with 304+ real SRA samples
- **Error scenarios**: Handles network failures, API errors, malformed data
- **Edge cases**: Manages missing fields, empty responses, timeouts

## 🔍 Troubleshooting

### Common Issues

1. **EDirect not found**: Ensure EDirect is installed and in PATH
2. **Proxy errors**: Script automatically unsets proxy variables
3. **Permission errors**: Ensure scripts have execute permissions
4. **Memory issues**: Process smaller batches for large datasets

### Debug Mode

Enable verbose logging:
```bash
export DEBUG=1
./build_metadata_enhanced.sh -i srr_list.txt -o output
```

## 📚 References

- [NCBI EDirect Documentation](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- [ENA API Documentation](https://ena-docs.readthedocs.io/en/latest/)
- [SRA Data Model](https://www.ncbi.nlm.nih.gov/sra/docs/sra_data_model/)
- [BioSample Documentation](https://www.ncbi.nlm.nih.gov/biosample/docs/)

## 🤝 Contributing

This pipeline is designed for production use with comprehensive error handling and documentation. Contributions are welcome for:
- Additional data source integrations
- Performance optimizations
- Enhanced error handling
- Documentation improvements

## 📄 License

This project is provided as-is for research and educational purposes. Please ensure compliance with NCBI and ENA data usage policies.

---

**Note**: This pipeline has been extensively tested and validated with real-world datasets. It represents a comprehensive solution for SRA metadata collection that goes beyond basic tools to provide complete, authoritative metadata from all available sources.