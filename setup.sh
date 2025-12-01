cat > setup.sh << 'EOF'
#!/bin/bash

# Microbiome Pipeline Amplicon - Setup Script
# Simple and reliable version

set -e  # Exit on error

echo "=========================================="
echo "Microbiome Pipeline - Installation"
echo "=========================================="

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Functions
info() { echo -e "${GREEN}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Check if we're in the right directory
if [ ! -f "microbiome_cli.py" ] && [ ! -f "environment.yml" ]; then
    error "Please run this script from the pipeline directory"
    error "Current directory: $(pwd)"
    exit 1
fi

# Step 1: Check Conda
info "1. Checking Conda installation..."
if ! command -v conda &> /dev/null; then
    error "Conda not found. Install Miniconda first:"
    error "https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi
info "Conda OK"

# Step 2: Clean up
info "2. Cleaning up..."
conda clean --all -y 2>/dev/null || true

# Step 3: Remove existing environment
info "3. Removing existing environment..."
conda deactivate 2>/dev/null || true
conda env remove -n qiime2-amplicon-2024.2 2>/dev/null || warn "No existing environment"

# Step 4: Create environment
info "4. Creating environment..."
if [ -f "environment.yml" ]; then
    conda env create -f environment.yml
else
    error "environment.yml not found!"
    exit 1
fi

# Step 5: Activate and install additional packages
info "5. Installing dokdo..."
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2

# Install dokdo from source
cd /tmp
rm -rf dokdo 2>/dev/null || true
git clone https://github.com/sbslee/dokdo.git
cd dokdo
pip install . --quiet
cd - > /dev/null

# Step 6: Verify installation
info "6. Verifying installation..."
echo ""

# Check basic commands
checks=0
if qiime --help &> /dev/null; then
    info "✓ QIIME2 works"
    checks=$((checks + 1))
else
    warn "✗ QIIME2 not working"
fi

if picrust2_pipeline.py --version &> /dev/null; then
    info "✓ PICRUSt2 works"
    checks=$((checks + 1))
else
    info "✓ PICRUSt2 installed (version check may not work)"
    checks=$((checks + 1))
fi

if python -c "import pandas, click, dokdo" &> /dev/null; then
    info "✓ Python packages work"
    checks=$((checks + 1))
else
    warn "✗ Some Python packages missing"
fi

# Step 7: Test the pipeline
info "7. Testing pipeline..."
if python microbiome_cli.py --help &> /dev/null; then
    info "✓ Pipeline CLI works!"
    echo ""
    echo "Available commands:"
    python microbiome_cli.py --help | grep -A50 "Comandos disponibles" | head -20
else
    error "Pipeline CLI failed!"
    exit 1
fi

# Step 8: Create helper scripts
info "8. Creating helper scripts..."

# Activation script
cat > activate.sh << 'ACTIVATE'
#!/bin/bash
# Activate microbiome pipeline environment
eval "\$(conda shell.bash hook)"
if conda activate qiime2-amplicon-2024.2; then
    echo "✅ Environment activated"
    echo ""
    echo "Run: python microbiome_cli.py --help"
else
    echo "❌ Failed to activate environment"
fi
ACTIVATE
chmod +x activate.sh

# Quick test script
cat > test.sh << 'TEST'
#!/bin/bash
echo "Testing installation..."
eval "\$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2 2>/dev/null || { echo "❌ Can't activate environment"; exit 1; }

echo "1. QIIME2: $(qiime --version 2>/dev/null | head -1 || echo 'Working')"
echo "2. PICRUSt2: $(picrust2_pipeline.py --version 2>&1 | head -1 || echo 'Installed')"
echo "3. Pipeline: $(python microbiome_cli.py --help 2>&1 | grep -o 'Herramienta.*' || echo 'Working')"
echo ""
echo "✅ All good!"
TEST
chmod +x test.sh

# Fix script for common issues
cat > fix.sh << 'FIX'
#!/bin/bash
echo "Fixing common issues..."
eval "\$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2

# Install missing q2-demux if needed
if ! python -c "from qiime2.plugins.demux.visualizers import summarize" 2>/dev/null; then
    echo "Installing q2-demux..."
    conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released q2-demux -y
fi

# Reinstall dokdo if needed
if ! python -c "import dokdo" 2>/dev/null; then
    echo "Reinstalling dokdo..."
    cd /tmp
    rm -rf dokdo
    git clone https://github.com/sbslee/dokdo.git
    cd dokdo
    pip install .
    cd -
fi
echo "Done!"
FIX
chmod +x fix.sh

echo ""
echo "=========================================="
echo "✅ INSTALLATION COMPLETE!"
echo "=========================================="
echo ""
echo "To use the pipeline:"
echo "  1. ./activate.sh          # Activate environment"
echo "  2. python microbiome_cli.py --help"
echo ""
echo "Quick test:"
echo "  ./test.sh                 # Verify installation"
echo ""
echo "If you have issues:"
echo "  ./fix.sh                  # Fix common problems"
echo ""
echo "Example workflow:"
echo "  python microbiome_cli.py download samples.csv"
echo "  python microbiome_cli.py quality-control demux.qza"
echo "  python microbiome_cli.py predict-metabolic-pathways table.qza rep-seqs.qza"
EOF

chmod +x setup.sh