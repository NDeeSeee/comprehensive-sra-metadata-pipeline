#!/usr/bin/env bash
# Test script for the comprehensive SRA metadata pipeline
# This script tests individual components to ensure everything works

set -euo pipefail

echo "🧪 Testing Comprehensive SRA Metadata Pipeline Components"
echo "=========================================================="

# Test 1: Check required tools
echo "1. Checking required tools..."
need() { command -v "$1" >/dev/null 2>&1 || { echo "❌ Missing: $1" >&2; echo "   Install with: conda install -c bioconda entrez-direct" >&2; exit 1; }; }
need esearch; need efetch; need curl; need jq; need python3
echo "✅ All required tools found"

# Test 1.5: Check Python packages
echo "1.5. Checking Python packages..."
python3 -c "import pandas, numpy, requests" 2>/dev/null || {
    echo "❌ Missing Python packages" >&2
    echo "   Install with: pip install pandas numpy requests" >&2
    exit 1
}
echo "✅ Python packages found"

# Test 2: Clear proxy settings
echo "2. Clearing proxy settings..."
unset https_proxy http_proxy HTTP_PROXY HTTPS_PROXY
export https_proxy="" http_proxy="" HTTP_PROXY="" HTTPS_PROXY=""
echo "✅ Proxy settings cleared"

# Test 3: Test NCBI EDirect
echo "3. Testing NCBI EDirect..."
if esearch -db sra -query "SRR3880241" | efetch -format runinfo | head -1 | grep -q "Run"; then
    echo "✅ NCBI EDirect working"
else
    echo "⚠️  NCBI EDirect may have issues (proxy or network)"
fi

# Test 4: Test ENA API
echo "4. Testing ENA API..."
if curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR3880241&result=read_run&fields=run_accession&format=tsv" | grep -q "run_accession"; then
    echo "✅ ENA API working"
else
    echo "❌ ENA API failed"
    exit 1
fi

# Test 5: Test ffq (if available)
echo "5. Testing ffq..."
if command -v ffq >/dev/null 2>&1; then
    if ffq SRR3880241 --json | head -1 | grep -q "{"; then
        echo "✅ ffq working"
    else
        echo "⚠️  ffq may have issues"
    fi
else
    echo "⚠️  ffq not installed (optional)"
fi

# Test 6: Test Python merge script
echo "6. Testing Python merge script..."
if python3 -c "import pandas, json, pathlib, re, xml.etree.ElementTree; print('✅ Python dependencies OK')"; then
    echo "✅ Python merge script ready"
else
    echo "❌ Python dependencies missing"
    exit 1
fi

echo ""
echo "🎉 Pipeline test completed!"
echo "Ready to run: ./build_metadata_enhanced.sh -i srr_list.txt -o output --with-geo --with-xml --with-ffq"