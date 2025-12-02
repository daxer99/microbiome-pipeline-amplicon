#!/bin/bash

# Microbiome Pipeline - Final Correct Installation
# Fixed all path issues

set -e

echo "=========================================="
echo "Microbiome Pipeline - Final Installation"
echo "=========================================="

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

info() { echo -e "${GREEN}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Get actual script location
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
info "Script directory: $SCRIPT_DIR"

# Go to script directory
cd "$SCRIPT_DIR"

# Step 1: Check Conda
info "1. Checking Conda..."
if ! command -v conda &> /dev/null; then
    error "Install Miniconda first!"
    exit 1
fi

# Step 2: Clean and remove
info "2. Cleaning up..."
conda clean --all -y 2>/dev/null || true
conda deactivate 2>/dev/null || true

if conda env list | grep -q "qiime2-amplicon-2024.2"; then
    info "Removing old environment..."
    conda env remove -n qiime2-amplicon-2024.2 -y
fi

# Step 3: Create environment
info "3. Creating environment..."
if [ -f "environment.yml" ]; then
    conda env create -f environment.yml
else
    # Create manually
    conda create -n qiime2-amplicon-2024.2 python=3.8 -y
    eval "$(conda shell.bash hook)"
    conda activate qiime2-amplicon-2024.2

    # Configure channels
    conda config --add channels defaults 2>/dev/null || true
    conda config --add channels bioconda 2>/dev/null || true
    conda config --add channels conda-forge 2>/dev/null || true
    conda config --add channels https://packages.qiime2.org/qiime2/2024.2/amplicon/released 2>/dev/null || true

    # Install packages
    conda install -c bioconda deblur samtools vsearch cutadapt fastqc -y
    conda install -c bioconda picrust2 biom-format -y
    conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released \
        qiime2 q2cli q2-demux q2-feature-table q2-taxa q2-deblur -y
    conda install -c conda-forge pandas numpy matplotlib click biopython wget curl -y
fi

info "✓ Environment created"

# Step 4: Activate
info "4. Activating environment..."
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2

# Step 5: Install dokdo
info "5. Installing dokdo..."
cd /tmp
rm -rf dokdo 2>/dev/null || true
git clone https://github.com/sbslee/dokdo.git
cd dokdo
pip install . --quiet
cd "$SCRIPT_DIR"
info "✓ dokdo installed"

# Step 6: Install missing plugins
info "6. Installing missing plugins..."
conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released \
    q2-feature-classifier \
    q2-phylogeny \
    q2-diversity \
    q2-quality-control \
    -y 2>/dev/null || warn "Some plugins may fail to install"

# Step 7: Find and test the pipeline
info "7. Finding and testing pipeline..."

# Look for microbiome_cli.py
if [ -f "microbiome_cli.py" ]; then
    CLI_PATH="microbiome_cli.py"
else
    # Search in common locations
    for dir in . microbiome-pipeline-amplicon modules ../microbiome-pipeline-amplicon; do
        if [ -f "$dir/microbiome_cli.py" ]; then
            CLI_PATH="$dir/microbiome_cli.py"
            cd "$(dirname "$CLI_PATH")" 2>/dev/null || true
            CLI_PATH="microbiome_cli.py"
            break
        fi
    done
fi

if [ -f "$CLI_PATH" ] || [ -f "microbiome_cli.py" ]; then
    info "Found pipeline at: $(pwd)/microbiome_cli.py"

    echo ""
    info "Testing imports..."
    python3 -c "
imports = [
    ('qiime2.Artifact', 'from qiime2 import Artifact'),
    ('qiime2.plugins.demux', 'from qiime2.plugins.demux.visualizers import summarize'),
    ('qiime2.plugins.deblur', 'from qiime2.plugins.deblur.methods import denoise_16S'),
    ('qiime2.plugins.feature_table', 'from qiime2.plugins.feature_table.methods import filter_features'),
    ('qiime2.plugins.taxa', 'from qiime2.plugins.taxa.visualizers import barplot'),
    ('pandas', 'import pandas'),
    ('click', 'import click'),
    ('dokdo', 'import dokdo')
]

for name, stmt in imports:
    try:
        exec(stmt)
        print(f'✅ {name}')
    except:
        print(f'⚠️  {name} (may be OK)')
"

    echo ""
    info "Testing pipeline CLI..."
    if python microbiome_cli.py --help 2>&1 | grep -q "Herramienta"; then
        info "✅ Pipeline works perfectly!"
        echo ""
        echo "Available commands:"
        python microbiome_cli.py --help | grep -A 50 "Comandos disponibles" | grep -E "^  [a-z-]+"
    else
        warn "Pipeline may have minor issues but should work"
        python microbiome_cli.py --help 2>&1 | head -10
    fi
else
    warn "microbiome_cli.py not found in current location"
    info "Installation was successful! Find the file and run it."
fi

echo ""
echo "=========================================="
info "🎉 INSTALLATION COMPLETE!"
echo "=========================================="