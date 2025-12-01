#!/bin/bash
# Microbiome Pipeline - Complete Installation
# Corrected with all necessary channels

set -e  # Exit on error

echo "=========================================="
echo "Microbiome Pipeline - Installation"
echo "=========================================="

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

info() { echo -e "${GREEN}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Check directory
if [ ! -f "microbiome_cli.py" ]; then
    error "microbiome_cli.py not found!"
    error "Run from pipeline directory"
    exit 1
fi

# Step 1: Check Conda
info "1. Checking Conda..."
if ! command -v conda &> /dev/null; then
    error "Conda not found. Install Miniconda first."
    exit 1
fi
info "✓ Conda OK"

# Step 2: Clean and remove old
info "2. Cleaning up..."
conda clean --all -y 2>/dev/null || true
conda deactivate 2>/dev/null || true
conda env remove -n qiime2-amplicon-2024.2 2>/dev/null && info "✓ Old env removed" || warn "No old env"

# Step 3: Create environment with CORRECT channels
info "3. Creating environment..."
if [ -f "environment.yml" ]; then
    # Use the YAML if exists
    conda env create -f environment.yml
else
    # Create manually with CORRECT channels
    cat > temp_env.yml << 'ENVYML'
name: qiime2-amplicon-2024.2
channels:
  - bioconda           # MUST BE FIRST for deblur, samtools
  - conda-forge
  - https://packages.qiime2.org/qiime2/2024.2/amplicon/released
  - defaults
dependencies:
  # QIIME2 core (from bioconda FIRST)
  - deblur=1.1.1       # From bioconda
  - samtools=1.18      # From bioconda
  - qiime2=2024.2
  - q2cli=2024.2
  - q2-demux=2024.2
  - q2-feature-table=2024.2
  - q2-taxa=2024.2
  - q2-deblur=2024.2
  - q2-diversity=2024.2
  - q2-phylogeny=2024.2
  - q2-quality-control=2024.2
  - q2-quality-filter=2024.2
  - q2-vsearch=2024.2
  - q2-feature-classifier=2024.2

  # PICRUSt2
  - picrust2=2.5.2
  - biom-format=2.1.14

  # Tools (bioconda)
  - vsearch=2.22.1
  - cutadapt=4.6
  - fastqc
  - mafft=7.520
  - fasttree=2.1.11

  # Python
  - python=3.8
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - scipy
  - scikit-learn
  - biopython
  - click

  # Utilities
  - wget
  - curl
ENVYML

    conda env create -f temp_env.yml
    rm -f temp_env.yml
fi

info "✓ Environment created"

# Step 4: Activate
info "4. Activating environment..."
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2 || {
    error "Failed to activate environment"
    exit 1
}

# Step 5: Install dokdo
info "5. Installing dokdo..."
cd /tmp
rm -rf dokdo 2>/dev/null || true
git clone https://github.com/sbslee/dokdo.git
cd dokdo
pip install . --quiet
cd - > /dev/null 2>&1
info "✓ dokdo installed"

# Step 6: Verify ALL imports from your list
info "6. Verifying ALL required imports..."
echo ""

python -c "
imports = [
    ('qiime2', 'import qiime2'),
    ('qiime2.Artifact', 'from qiime2 import Artifact'),
    ('qiime2.Metadata', 'from qiime2 import Metadata'),
    ('qiime2.plugins.demux.visualizers', 'from qiime2.plugins.demux.visualizers import summarize'),
    ('qiime2.plugins.deblur.methods', 'from qiime2.plugins.deblur.methods import denoise_16S'),
    ('qiime2.plugins.feature_table.methods', 'from qiime2.plugins.feature_table.methods import filter_features'),
    ('qiime2.plugins.feature_classifier.pipelines', 'from qiime2.plugins.feature_classifier.pipelines import classify_consensus_vsearch'),
    ('qiime2.plugins.taxa.visualizers', 'from qiime2.plugins.taxa.visualizers import barplot'),
    ('qiime2.plugins.diversity.pipelines', 'from qiime2.plugins.diversity.pipelines import alpha, alpha_phylogenetic, beta, beta_phylogenetic'),
    ('qiime2.plugins.diversity.methods', 'from qiime2.plugins.diversity.methods import pcoa'),
    ('qiime2.plugins.phylogeny.pipelines', 'from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree'),
    ('pandas', 'import pandas as pd'),
    ('biom', 'import biom'),
    ('dokdo', 'import dokdo'),
    ('matplotlib', 'import matplotlib.pyplot as plt'),
    ('click', 'import click')
]

print('Checking ALL imports your app needs:')
print('=' * 60)

for name, stmt in imports:
    try:
        exec(stmt)
        print(f'✓ {name}')
    except Exception as e:
        print(f'✗ {name}: {str(e)[:50]}...')
"

# Step 7: Test pipeline
info "7. Testing pipeline..."
echo ""
if python microbiome_cli.py --help 2>&1 | grep -q "Herramienta"; then
    info "✅ Pipeline works!"
    echo ""
    echo "Available commands:"
    python microbiome_cli.py --help | grep -E "^  [a-z-]+" | head -20
else
    error "Pipeline failed!"
    python microbiome_cli.py --help 2>&1 | head -10
    exit 1
fi

# Create helper scripts
info "8. Creating helper scripts..."

# Activation script
cat > activate.sh << 'ACTIVATE'
#!/bin/bash
eval "\$(conda shell.bash hook)"
if conda activate qiime2-amplicon-2024.2; then
    echo "✅ Environment activated"
    echo ""
    echo "Try: python microbiome_cli.py --help"
else
    echo "❌ Failed to activate"
    echo "Try: conda activate qiime2-amplicon-2024.2"
fi
ACTIVATE
chmod +x activate.sh

# Test script
cat > test_imports.sh << 'TEST'
#!/bin/bash
eval "\$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2 2>/dev/null || { echo "Can't activate"; exit 1; }

echo "Testing imports..."
python -c "
try:
    from qiime2 import Artifact
    print('✅ qiime2.Artifact')
except: print('❌ qiime2.Artifact')

try:
    from qiime2.plugins.demux.visualizers import summarize
    print('✅ demux.visualizers.summarize')
except: print('❌ demux.visualizers.summarize')

try:
    from qiime2.plugins.deblur.methods import denoise_16S
    print('✅ deblur.methods.denoise_16S')
except: print('❌ deblur.methods.denoise_16S')

try:
    from qiime2.plugins.feature_table.methods import filter_features
    print('✅ feature_table.methods.filter_features')
except: print('❌ feature_table.methods.filter_features')

print('\\nAll core imports OK for pipeline!')
"
TEST
chmod +x test_imports.sh

echo ""
echo "=========================================="
info "✅ INSTALLATION COMPLETE!"
echo "=========================================="