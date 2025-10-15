# Installation Guide - Cancer Analysis Pipeline Dependencies

## 🎯 **Quick Installation**

The Cancer Analysis Pipeline requires two main dependencies to function completely:

### **1. Entrez Direct Tools**
Required for metadata collection from NCBI databases.

```bash
# Install via conda (recommended)
conda install -c bioconda entrez-direct

# Or install via apt (Ubuntu/Debian)
sudo apt-get install ncbi-entrez-direct

# Verify installation
esearch --version
efetch --version
```

### **2. SRA Toolkit**
Required for FASTQ file downloads.

```bash
# Install via conda (recommended)
conda install -c bioconda sra-tools

# Or download from NCBI
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
export PATH=$PATH:$(pwd)/sratoolkit.current-ubuntu64/bin

# Verify installation
prefetch --version
fastq-dump --version
```

## 🧪 **Test Installation**

After installation, test the complete pipeline:

```bash
# Test cancer type search (works without dependencies)
python3 scripts/cancer_type_search.py -c "lung cancer" --test

# Test complete pipeline (requires dependencies)
python3 scripts/cancer_analysis_pipeline.py -c "lung cancer" --test
```

## 📋 **System Requirements**

- **Python 3.8+** with pandas, requests, json
- **Unix/Linux environment** with bash, curl, jq
- **Internet connection** for database access
- **Sufficient disk space** for FASTQ downloads (GB to TB)

## 🔧 **Troubleshooting**

### **Common Issues**

1. **"Missing: esearch"**
   ```bash
   conda install -c bioconda entrez-direct
   ```

2. **"Missing: prefetch"**
   ```bash
   conda install -c bioconda sra-tools
   ```

3. **Permission errors**
   ```bash
   chmod +x scripts/*.py scripts/*.sh
   ```

4. **Python path issues**
   ```bash
   which python3
   # Use full path if needed: /usr/local/anaconda3-2020/bin/python3
   ```

### **Verification Commands**

```bash
# Check all dependencies
which python3 esearch efetch prefetch fastq-dump curl jq

# Test Python packages
python3 -c "import pandas, requests, json; print('All packages available')"

# Test pipeline components
python3 scripts/cancer_type_search.py --help
python3 scripts/cancer_classification.py --help
```

## 🚀 **Ready to Use**

Once dependencies are installed, the pipeline is ready for production use:

```bash
# Example: Complete lung cancer analysis
python3 scripts/cancer_analysis_pipeline.py -c "lung cancer" -o lung_cancer_results

# Example: Breast cancer with limited downloads
python3 scripts/cancer_analysis_pipeline.py -c "breast cancer" --test --max-per-category 5
```

The system will automatically:
1. Search for cancer type and discover SRR IDs
2. Collect comprehensive metadata from multiple sources
3. Classify samples into Tumor/Normal/Cell line categories
4. Download and organize FASTQ files by category
5. Generate summary reports and logs