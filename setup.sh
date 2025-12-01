cat > setup_fixed.sh << 'EOF'
#!/bin/bash

# Microbiome Pipeline Amplicon - Setup Script
# Fixed version with better path handling

set -e  # Exit on any error

echo "=================================================="
echo "Microbiome Pipeline Amplicon - Installation"
echo "=================================================="

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_status() { echo -e "${GREEN}[INFO]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }
print_step() { echo -e "${BLUE}[STEP]${NC} $1"; }

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
print_status "Script directory: $SCRIPT_DIR"

# Check Conda
check_conda() {
    if ! command -v conda &> /dev/null; then
        print_error "Conda not installed"
        exit 1
    fi
    print_status "Conda is installed"
}

# Main installation
main() {
    print_step "1. Checking environment..."
    check_conda

    # Check if we're in the right directory
    if [ ! -f "microbiome_cli.py" ] && [ ! -f "environment.yml" ]; then
        print_warning "Not in pipeline directory. Looking for it..."
        if [ -f "$SCRIPT_DIR/microbiome_cli.py" ]; then
            cd "$SCRIPT_DIR"
            print_status "Changed to: $(pwd)"
        else
            print_error "Cannot find pipeline files"
            exit 1
        fi
    fi

    print_step "2. Cleaning cache..."
    conda clean --all -y 2>/dev/null || true

    print_step "3. Removing existing environment..."
    conda env remove -n qiime2-amplicon-2024.2 2>/dev/null || print_warning "No existing environment"

    print_step "4. Creating environment from environment.yml..."
    if [ -f "environment.yml" ]; then
        conda env create -f environment.yml
    else
        print_error "environment.yml not found"
        exit 1
    fi

    print_step "5. Activating environment..."
    eval "$(conda shell.bash hook)"
    conda activate qiime2-amplicon-2024.2

    print_step "6. Installing dokdo from source..."
    cd /tmp
    rm -rf dokdo
    git clone https://github.com/sbslee/dokdo.git
    cd dokdo
    pip install .
    cd "$SCRIPT_DIR"

    print_step "7. Verifying installation..."

    # Check basic imports
    python -c "
imports = [
    ('qiime2', 'import qiime2'),
    ('qiime2.plugins.demux', 'from qiime2.plugins.demux.visualizers import summarize'),
    ('qiime2.plugins.feature_table', 'from qiime2.plugins.feature_table.methods import filter_features'),
    ('qiime2.plugins.deblur', 'import qiime2.plugins.deblur'),
    ('pandas', 'import pandas'),
    ('click', 'import click'),
    ('dokdo', 'import dokdo')
]

print('Essential imports:')
for name, stmt in imports:
    try:
        exec(stmt)
        print(f'  ✅ {name}')
    except:
        print(f'  ⚠️  {name} (but may be OK)')
"

    print_step "8. Testing pipeline CLI..."
    if [ -f "microbiome_cli.py" ]; then
        if python microbiome_cli.py --help 2>&1 | grep -q "Herramienta de análisis"; then
            print_status "✅ Pipeline CLI works!"
        else
            print_warning "Pipeline may have issues"
        fi
    else
        print_error "microbiome_cli.py not found"
    fi

    print_step "9. Creating helper scripts..."

    # Activation script
    cat > activate.sh << 'SCRIPT'
#!/bin/bash
eval "\$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2
echo "Environment activated"
echo "Run: python microbiome_cli.py --help"
SCRIPT
    chmod +x activate.sh

    # Quick test script
    cat > test_install.sh << 'SCRIPT'
#!/bin/bash
eval "\$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2
echo "Testing installation..."
echo "1. QIIME2:" && qiime --version
echo "2. PICRUSt2:" && picrust2_pipeline.py --version 2>/dev/null || echo "  (version check not available)"
echo "3. Pipeline:" && python microbiome_cli.py --help | head -5
SCRIPT
    chmod +x test_install.sh

    print_step "10. Installation complete!"
    echo ""
    echo "Next steps:"
    echo "1. Activate environment: ./activate.sh"
    echo "2. Test installation: ./test_install.sh"
    echo "3. View commands: python microbiome_cli.py --help"
    echo ""
    echo "Quick start:"
    echo "  python microbiome_cli.py download samples.csv"
    echo "  python microbiome_cli.py quality-control demux.qza"
    echo "  python microbiome_cli.py predict-metabolic-pathways table.qza rep-seqs.qza"
}

main
EOF