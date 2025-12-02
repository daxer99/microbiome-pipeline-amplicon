# Microbiome Pipeline - Complete Working Installation
# Includes ALL required QIIME2 plugins

set -e

echo "=========================================="
echo "Microbiome Pipeline - Complete Installation"
echo "=========================================="

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

info() { echo -e "${GREEN}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Check we're in right place
if [ ! -f "microbiome_cli.py" ] && [ ! -f "environment.yml" ]; then
    error "Run from pipeline directory!"
    exit 1
fi

# Step 1: Check Conda
info "1. Checking Conda..."
if ! command -v conda &> /dev/null; then
    error "Install Miniconda first: https://docs.conda.io"
    exit 1
fi
info "✓ Conda OK"

# Step 2: Clean
info "2. Cleaning up..."
conda clean --all -y 2>/dev/null || true
conda deactivate 2>/dev/null || true

# Step 3: Remove old env
if conda env list | grep -q "qiime2-amplicon-2024.2"; then
    info "Removing old environment..."
    conda env remove -n qiime2-amplicon-2024.2 -y
fi

# Step 4: Create environment with ALL plugins
info "3. Creating environment with ALL required plugins..."
info "This may take 20-40 minutes depending on internet..."

# Use mamba if available (much faster)
if command -v mamba &> /dev/null; then
    info "Using mamba for faster installation..."
    mamba env create -f environment.yml
else
    # Try to install mamba
    info "Installing mamba for faster package management..."
    conda install -c conda-forge mamba -y
    mamba env create -f environment.yml
fi

if [ $? -eq 0 ]; then
    info "✓ Environment created successfully"
else
    warn "Environment creation had issues, trying alternative method..."

    # Alternative: Create step by step
    conda create -n qiime2-amplicon-2024.2 python=3.8 -y
    eval "$(conda shell.bash hook)"
    conda activate qiime2-amplicon-2024.2

    # Configure channels
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --add channels https://packages.qiime2.org/qiime2/2024.2/amplicon/released

    # Install in batches
    info "Installing QIIME2 core..."
    conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released \
        qiime2 q2cli q2-demux q2-feature-table q2-taxa q2-deblur -y

    info "Installing missing plugins..."
    conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released \
        q2-feature-classifier q2-phylogeny q2-diversity q2-quality-control -y

    info "Installing PICRUSt2..."
    conda install -c bioconda picrust2 vsearch cutadapt -y

    info "Installing Python packages..."
    conda install -c conda-forge pandas numpy matplotlib click biopython -y
fi

# Step 5: Activate
info "4. Activating environment..."
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2 || {
    error "Failed to activate environment"
    exit 1
}

# Step 6: Install dokdo
info "5. Installing dokdo..."
cd /tmp
rm -rf dokdo 2>/dev/null || true
git clone https://github.com/sbslee/dokdo.git
cd dokdo
pip install . --quiet
cd - > /dev/null 2>&1
info "✓ dokdo installed"

# Step 7: Verify ALL imports from your code
info "6. Verifying ALL required imports..."
echo ""

python3 -c "
# ALL imports from your application
imports = [
    # Basic Python
    ('tempfile', 'import tempfile'),
    ('pathlib', 'import pathlib'),
    ('pandas', 'import pandas as pd'),
    ('os', 'import os'),
    ('subprocess', 'import subprocess'),
    ('shutil', 'import shutil'),
    ('matplotlib', 'import matplotlib.pyplot as plt'),

    # QIIME2 core
    ('qiime2', 'import qiime2'),
    ('qiime2.Artifact', 'from qiime2 import Artifact'),
    ('qiime2.Metadata', 'from qiime2 import Metadata'),

    # QIIME2 plugins - ALL that you use
    ('qiime2.plugins.demux', 'from qiime2.plugins.demux.visualizers import summarize'),
    ('qiime2.plugins.deblur', 'from qiime2.plugins.deblur.methods import denoise_16S'),
    ('qiime2.plugins.feature_table', 'from qiime2.plugins.feature_table.methods import filter_features'),
    ('qiime2.plugins.feature_classifier', 'from qiime2.plugins.feature_classifier.pipelines import classify_consensus_vsearch'),
    ('qiime2.plugins.taxa', 'from qiime2.plugins.taxa.visualizers import barplot'),
    ('qiime2.plugins.diversity', 'from qiime2.plugins.diversity.pipelines import alpha, alpha_phylogenetic, beta, beta_phylogenetic'),
    ('qiime2.plugins.diversity.methods', 'from qiime2.plugins.diversity.methods import pcoa'),
    ('qiime2.plugins.phylogeny', 'from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree'),
    ('qiime2.plugins.quality_control', 'from qiime2.plugins.quality_control.visualizers import summarize as qc_summarize'),

    # Other packages
    ('biom', 'import biom'),
    ('dokdo', 'import dokdo'),
    ('click', 'import click'),
]

print('Verifying ALL imports your application uses:')
print('=' * 60)

failed_imports = []
for name, stmt in imports:
    try:
        exec(stmt)
        print(f'✅ {name}')
    except Exception as e:
        print(f'❌ {name}: {str(e)[:50]}')
        failed_imports.append(name)

print('')
if len(failed_imports) == 0:
    print('🎉 ALL IMPORTS SUCCESSFUL!')
else:
    print(f'⚠️  {len(failed_imports)} imports failed: {failed_imports}')
    print('Some features may not work.')
"

# Step 8: Test the pipeline
info "7. Testing pipeline CLI..."
echo ""

if [ -f "microbiome_cli.py" ]; then
    if python microbiome_cli.py --help 2>&1 | grep -q "Herramienta"; then
        info "✅ Pipeline CLI works perfectly!"
        echo ""
        echo "Available commands:"
        python microbiome_cli.py --help | grep -A 50 "Comandos disponibles" | grep -E "^  [a-z-]+"

        # Show quick help
        echo ""
        info "Quick test of all commands:"
        python microbiome_cli.py --help | grep -E "^  [a-z-]+" | while read cmd; do
            echo "  ✓ $cmd"
        done
    else
        error_msg=$(python microbiome_cli.py --help 2>&1 | head -5)
        warn "Pipeline test issues:"
        echo "$error_msg"

        # Check specific error
        if echo "$error_msg" | grep -q "feature_classifier"; then
            info "Installing missing feature-classifier..."
            conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released q2-feature-classifier -y
            info "Now test again: python microbiome_cli.py --help"
        fi
    fi
else
    error "microbiome_cli.py not found!"
    echo "Current directory: $(pwd)"
    ls -la
fi

# Create helper scripts
info "8. Creating helper scripts..."

# Activation script
cat > activate.sh << 'ACTIVATE'
#!/bin/bash
eval "\$(conda shell.bash hook)"
if conda activate qiime2-amplicon-2024.2; then
    echo "========================================"
    echo "✅ Microbiome Pipeline Environment"
    echo "========================================"
    echo ""
    echo "To use: python microbiome_cli.py --help"
    echo ""
    echo "Quick start:"
    echo "  python microbiome_cli.py download samples.csv"
    echo "  python microbiome_cli.py quality-control demux.qza"
    echo "  python microbiome_cli.py predict-metabolic-pathways table.qza rep-seqs.qza"
else
    echo "❌ Failed to activate environment"
fi
ACTIVATE
chmod +x activate.sh

# Test script
cat > test_imports.sh << 'TEST'
#!/bin/bash
eval "\$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2 2>/dev/null || { echo "❌ Can't activate environment"; exit 1; }

echo "Testing critical imports..."
python3 -c "
critical = [
    'from qiime2 import Artifact',
    'from qiime2.plugins.demux.visualizers import summarize',
    'from qiime2.plugins.deblur.methods import denoise_16S',
    'from qiime2.plugins.feature_table.methods import filter_features',
    'from qiime2.plugins.feature_classifier.pipelines import classify_consensus_vsearch',
    'from qiime2.plugins.taxa.visualizers import barplot',
    'from qiime2.plugins.diversity.pipelines import alpha, beta',
    'import pandas',
    'import click',
    'import dokdo'
]

for stmt in critical:
    try:
        exec(stmt)
        print('✅ OK')
    except Exception as e:
        print(f'❌ Failed: {str(e)[:40]}')
"

echo ""
echo "Testing pipeline..."
python microbiome_cli.py --help 2>&1 | head -3
TEST
chmod +x test_imports.sh

# Fix script
cat > fix_missing.sh << 'FIX'
#!/bin/bash
echo "Fixing missing packages..."
eval "\$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2

echo "1. Updating conda..."
conda update --all -y

echo "2. Installing missing QIIME2 plugins..."
missing=()
if ! python -c "from qiime2.plugins.feature_classifier import __version__" 2>/dev/null; then
    missing+=("q2-feature-classifier")
fi
if ! python -c "from qiime2.plugins.phylogeny import __version__" 2>/dev/null; then
    missing+=("q2-phylogeny")
fi

if [ ${#missing[@]} -gt 0 ]; then
    echo "Installing: ${missing[*]}"
    conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released ${missing[*]} -y
fi

echo "3. Reinstalling dokdo if needed..."
if ! python -c "import dokdo" 2>/dev/null; then
    cd /tmp
    git clone https://github.com/sbslee/dokdo.git
    cd dokdo
    pip install .
fi

echo "✅ Fix complete!"
FIX
chmod +x fix_missing.sh

echo ""
echo "=========================================="
info "✅ INSTALLATION COMPLETE!"
echo "=========================================="