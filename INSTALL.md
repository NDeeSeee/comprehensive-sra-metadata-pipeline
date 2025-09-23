# Installation Guide for Comprehensive SRA Metadata Pipeline

## 🚀 Complete Installation Instructions

### Option 1: Conda Environment (Recommended)

Create a dedicated conda environment with all required tools:

```bash
# Create conda environment
conda create -n sra-metadata python=3.9 -y
conda activate sra-metadata

# Install Python dependencies
conda install pandas numpy requests -y

# Install NCBI EDirect tools
conda install -c bioconda entrez-direct -y

# Install ffq (optional but recommended)
pip install ffq

# Install additional tools
conda install curl jq -y
```

### Option 2: Manual Installation

If you prefer to install tools system-wide or use existing installations:

#### 1. Install NCBI EDirect Tools
```bash
# Download and install EDirect
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# Add to PATH (add to ~/.bashrc or ~/.zshrc)
echo 'export PATH=$PATH:$HOME/edirect' >> ~/.bashrc
source ~/.bashrc
```

#### 2. Install Python Dependencies
```bash
pip install pandas numpy requests
```

#### 3. Install ffq (Optional)
```bash
pip install ffq
```

#### 4. Install System Tools
```bash
# On Ubuntu/Debian
sudo apt-get install curl jq

# On macOS
brew install curl jq

# On CentOS/RHEL
sudo yum install curl jq
```

### Option 3: Docker (Alternative)

If you prefer containerized deployment:

```dockerfile
FROM ubuntu:20.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    jq \
    python3 \
    python3-pip \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip3 install pandas numpy requests ffq

# Install EDirect
RUN sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# Set PATH
ENV PATH="/root/edirect:${PATH}"

# Set working directory
WORKDIR /app

# Copy pipeline files
COPY . .

# Make scripts executable
RUN chmod +x scripts/*.sh
```

## 🔧 Post-Installation Setup

### 1. Configure NCBI API Key (Optional but Recommended)
```bash
# Get API key from: https://www.ncbi.nlm.nih.gov/account/settings/
export NCBI_API_KEY="your_api_key_here"

# Add to ~/.bashrc for persistence
echo 'export NCBI_API_KEY="your_api_key_here"' >> ~/.bashrc
```

### 2. Test Installation
```bash
# Clone the repository
git clone https://github.com/NDeeSeee/comprehensive-sra-metadata-pipeline.git
cd comprehensive-sra-metadata-pipeline

# Test all components
cd scripts
./test_pipeline.sh
```

### 3. Verify Individual Tools
```bash
# Test EDirect
esearch -db sra -query "SRR1234567" | efetch -format runinfo | head -1

# Test ENA API
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR1234567&result=read_run&fields=run_accession&format=tsv"

# Test ffq
ffq SRR1234567 --json | head -1

# Test Python
python3 -c "import pandas, numpy, requests; print('Python dependencies OK')"
```

## 🎯 Quick Start Commands

### Using Conda (Recommended)
```bash
# Create and activate environment
conda create -n sra-metadata python=3.9 -y
conda activate sra-metadata

# Install all dependencies
conda install -c bioconda entrez-direct pandas numpy requests curl jq -y
pip install ffq

# Clone and run
git clone https://github.com/NDeeSeee/comprehensive-sra-metadata-pipeline.git
cd comprehensive-sra-metadata-pipeline/scripts
./test_pipeline.sh
```

### Using System Installation
```bash
# Install EDirect
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
echo 'export PATH=$PATH:$HOME/edirect' >> ~/.bashrc
source ~/.bashrc

# Install Python dependencies
pip install pandas numpy requests ffq

# Install system tools (Ubuntu/Debian example)
sudo apt-get install curl jq

# Clone and run
git clone https://github.com/NDeeSeee/comprehensive-sra-metadata-pipeline.git
cd comprehensive-sra-metadata-pipeline/scripts
./test_pipeline.sh
```

## 🔍 Troubleshooting

### Common Issues

1. **EDirect not found**: Ensure EDirect is installed and in PATH
2. **Permission denied**: Make scripts executable with `chmod +x scripts/*.sh`
3. **Python module not found**: Install missing modules with `pip install module_name`
4. **curl/jq not found**: Install system tools for your OS

### Environment Issues

If you have multiple Python environments:
```bash
# Check which Python is being used
which python3
python3 --version

# Ensure you're using the right environment
conda activate sra-metadata  # if using conda
```

### Proxy Issues (if any)
```bash
# Clear proxy settings
unset https_proxy http_proxy HTTP_PROXY HTTPS_PROXY
export https_proxy="" http_proxy="" HTTP_PROXY="" HTTPS_PROXY=""
```

## 📋 Requirements Summary

**Required Tools:**
- `esearch`, `efetch` (NCBI EDirect)
- `curl` (system tool)
- `jq` (system tool)
- `python3` (Python 3.7+)
- `pandas`, `numpy`, `requests` (Python packages)

**Optional Tools:**
- `ffq` (additional metadata source)
- NCBI API key (for higher rate limits)

**Minimum System Requirements:**
- 4GB RAM (for processing large datasets)
- 1GB disk space (for output files)
- Internet connection (for API calls)