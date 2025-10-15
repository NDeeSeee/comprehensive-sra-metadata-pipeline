# Setup Instructions

## Quick Setup

1. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Install Entrez Direct:**
   ```bash
   # macOS
   brew install ncbi-entrez-direct
   
   # Linux
   wget https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
   tar -xzf edirect.tar.gz
   export PATH=$PATH:$PWD/edirect
   ```

3. **Run the pipeline:**
   ```bash
   ./scripts/comprehensive_metadata_pipeline.sh -i data/srr_list.txt --classify
   ```

## Optional Tools

- **ffq** (for file information): `pip install ffq`
- **jq** (for JSON processing): `brew install jq` or `apt-get install jq`
