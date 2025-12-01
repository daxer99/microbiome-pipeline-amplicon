#!/bin/bash

# Microbiome Pipeline Amplicon - Setup Script
# Updated with q2-demux and proper dokdo installation

set -e  # Exit on any error

echo "=================================================="
echo "Microbiome Pipeline Amplicon - Complete Installation"
echo "=================================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

# Check if Conda is installed
check_conda() {
    if ! command -v conda &> /dev/null; then
        print_error "Conda is not installed. Please install Miniconda or Anaconda first."
        echo "Download Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
    print_status "Conda is installed"
}

# Clean Conda cache
clean_conda_cache() {
    print_step "Cleaning Conda cache..."
    conda clean --all -y
}

# Remove existing environment
remove_existing_env() {
    local env_name="qiime2-amplicon-2024.2"

    if conda env list | grep -q "$env_name"; then
        print_warning "Environment $env_name already exists"
        read -p "Do you want to remove and recreate it? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            conda env remove -n "$env_name"
            print_status "Environment removed"
        else
            print_status "Using existing environment"
            return 1  # Don't recreate
        fi
    fi
    return 0  # Recreate
}

# Create Conda environment
create_environment() {
    local env_name="qiime2-amplicon-2024.2"

    print_step "Creating Conda environment: $env_name"

    # Create from YAML if exists
    if [ -f "environment.yml" ]; then
        print_status "Using environment.yml file..."
        conda env create -f environment.yml
    else
        print_status "Creating environment manually..."

        # Create base environment
        conda create -n "$env_name" python=3.8 -y

        # Activate
        eval "$(conda shell.bash hook)"
        conda activate "$env_name"

        # Install QIIME2 with ALL essential plugins (including q2-demux)
        print_status "Installing QIIME2 with essential plugins..."
        conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released \
            qiime2 \
            q2cli \
            q2-demux \           # <-- CRITICAL for quality_control.py
            q2-feature-table \
            q2-taxa \
            q2-deblur \
            q2-dada2 \
            q2-quality-control \
            q2-quality-filter \
            q2-vsearch \
            q2-feature-classifier \
            q2-diversity \
            q2-phylogeny \
            -y

        # Install PICRUSt2
        print_status "Installing PICRUSt2..."
        conda install -c bioconda picrust2 -y

        # Install bioinformatics tools
        print_status "Installing bioinformatics tools..."
        conda install -c bioconda vsearch cutadapt fastqc multiqc -y

        # Install Python packages
        print_status "Installing Python packages..."
        conda install -c conda-forge pandas numpy matplotlib seaborn scipy scikit-learn click biopython -y

        # Install utilities
        print_status "Installing utilities..."
        conda install -c conda-forge wget curl parallel pigz pbzip2 -y
    fi

    print_status "Environment created successfully"
}

# Install dokdo from source (as per author's instructions)
install_dokdo() {
    print_step "Installing dokdo from source..."

    # Check if already installed
    if python -c "import dokdo" &> /dev/null; then
        print_status "dokdo is already installed"
        return 0
    fi

    # Clone and install from source
    local temp_dir=$(mktemp -d)
    cd "$temp_dir"

    print_status "Cloning dokdo repository..."
    git clone https://github.com/sbslee/dokdo.git
    cd dokdo

    print_status "Installing dokdo..."
    pip install .

    # Clean up
    cd ~
    rm -rf "$temp_dir"

    # Verify installation
    if python -c "import dokdo" &> /dev/null; then
        print_status "dokdo installed successfully"
    else
        print_warning "dokdo installation may have failed"
        print_status "You can install it manually: git clone https://github.com/sbslee/dokdo && cd dokdo && pip install ."
    fi
}

# Verify QIIME2 Python imports
verify_qiime2_imports() {
    print_step "Verifying QIIME2 Python imports..."

    python -c "
imports_to_check = [
    ('qiime2', 'import qiime2'),
    ('qiime2.plugins.demux', 'from qiime2.plugins.demux.visualizers import summarize'),
    ('qiime2.plugins.feature_table', 'from qiime2.plugins.feature_table.methods import filter_features'),
    ('qiime2.plugins.quality_control', 'from qiime2.plugins.quality_control.visualizers import summarize'),
    ('qiime2.plugins.quality_filter', 'from qiime2.plugins.quality_filter.methods import q_score'),
    ('qiime2.plugins.deblur', 'from qiime2.plugins.deblur.methods import denoise_16S'),
    ('qiime2.plugins.taxa', 'from qiime2.plugins.taxa.methods import classify_sklearn'),
]

print('Checking critical imports for the pipeline:')
print('=' * 50)

all_ok = True
for name, import_stmt in imports_to_check:
    try:
        exec(import_stmt)
        print(f'✅ {name}')
    except ImportError as e:
        print(f'❌ {name}: {e}')
        all_ok = False

if all_ok:
    print('\\n✅ All critical imports are working!')
else:
    print('\\n⚠️  Some imports failed. The pipeline may not work correctly.')
"
}

# Test the pipeline
test_pipeline() {
    print_step "Testing the pipeline..."

    if python microbiome_cli.py --help 2>&1 | grep -q "CLI Principal"; then
        print_status "✅ Pipeline CLI is working!"
        echo ""
        python microbiome_cli.py --help | head -30
    else
        print_error "Pipeline CLI is not working correctly"
        print_status "Checking for errors..."
        python microbiome_cli.py --help 2>&1 | head -20
    fi
}

# Create activation script
create_activation_script() {
    cat > activate_pipeline.sh << 'EOF'
#!/bin/bash
# Script to activate the microbiome pipeline environment

echo "=================================================="
echo "Microbiome Pipeline Environment Activation"
echo "=================================================="

# Initialize conda
eval "$(conda shell.bash hook)"

# Activate environment
if conda activate qiime2-amplicon-2024.2; then
    echo ""
    echo "✅ Environment activated successfully!"
    echo ""
    echo "Available tools:"
    echo "  • QIIME2: qiime --help"
    echo "  • PICRUSt2: picrust2_pipeline.py --version"
    echo "  • vsearch: vsearch --version"
    echo "  • cutadapt: cutadapt --version"
    echo ""
    echo "Pipeline commands:"
    echo "  python microbiome_cli.py --help"
    echo ""
    echo "Quick start:"
    echo "  1. Download data: python microbiome_cli.py download samples.csv"
    echo "  2. Quality control: python microbiome_cli.py quality-control demux.qza"
    echo "  3. Run analysis: python microbiome_cli.py predict-metabolic-pathways table.qza rep-seqs.qza"
else
    echo ""
    echo "❌ Failed to activate environment."
    echo "Try: conda activate qiime2-amplicon-2024.2"
fi
EOF

    chmod +x activate_pipeline.sh
    print_status "Created activation script: ./activate_pipeline.sh"
}

# Create verification script
create_verification_script() {
    cat > verify_installation.sh << 'EOF'
#!/bin/bash
echo "Verifying Microbiome Pipeline Installation"
echo "=========================================="

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2 2>/dev/null

if [ $? -ne 0 ]; then
    echo "❌ Cannot activate environment"
    exit 1
fi

echo ""
echo "1. Checking QIIME2 plugins..."
qiime info | grep -A10 "Installed plugins"

echo ""
echo "2. Checking critical tools..."
tools=("qiime" "picrust2_pipeline.py" "vsearch" "cutadapt" "fastqc")
for tool in "${tools[@]}"; do
    if command -v $tool &> /dev/null; then
        echo "  ✅ $tool"
    else
        echo "  ❌ $tool (MISSING)"
    fi
done

echo ""
echo "3. Checking Python imports..."
python -c "
imports = [
    'qiime2',
    'qiime2.plugins.demux',
    'qiime2.plugins.feature_table',
    'qiime2.plugins.quality_control',
    'pandas',
    'click',
    'biopython'
]

for module in imports:
    try:
        __import__(module)
        print(f'  ✅ {module}')
    except ImportError as e:
        print(f'  ❌ {module}: {e}')
"

echo ""
echo "4. Testing pipeline..."
if python microbiome_cli.py --help &> /dev/null; then
    echo "  ✅ Pipeline CLI works"
else
    echo "  ❌ Pipeline CLI failed"
fi

echo ""
echo "=========================================="
echo "Verification complete!"
EOF

    chmod +x verify_installation.sh
    print_status "Created verification script: ./verify_installation.sh"
}

# Create troubleshooting script
create_troubleshooting_script() {
    cat > fix_common_issues.sh << 'EOF'
#!/bin/bash
echo "Fixing Common Installation Issues"
echo "================================="

# Activate environment
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2

echo ""
echo "1. Checking for missing QIIME2 plugins..."
missing_plugins=()

# Check each critical plugin
python -c "
plugins = ['demux', 'feature_table', 'quality_control', 'quality_filter', 'deblur', 'taxa']
import qiime2.plugins as q2p

for plugin in plugins:
    try:
        getattr(q2p, plugin)
        print(f'✅ q2-{plugin} is installed')
    except AttributeError:
        print(f'❌ q2-{plugin} is MISSING')
        missing_plugins.append(f'q2-{plugin}')
"

echo ""
echo "2. Installing missing plugins..."
if [ ${#missing_plugins[@]} -gt 0 ]; then
    echo "Installing: ${missing_plugins[*]}"
    conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released ${missing_plugins[*]} -y
else
    echo "No missing plugins found"
fi

echo ""
echo "3. Checking dokdo..."
if python -c "import dokdo" &> /dev/null; then
    echo "✅ dokdo is installed"
else
    echo "Installing dokdo from source..."
    cd /tmp
    git clone https://github.com/sbslee/dokdo.git
    cd dokdo
    pip install .
    echo "✅ dokdo installed"
fi

echo ""
echo "4. Final test..."
python microbiome_cli.py --help | head -5
EOF

    chmod +x fix_common_issues.sh
    print_status "Created troubleshooting script: ./fix_common_issues.sh"
}

# Main installation function
main() {
    echo
    print_step "Starting complete installation process..."

    # Check conda
    check_conda

    # Clean cache
    clean_conda_cache

    # Remove existing env if needed
    if remove_existing_env; then
        # Create environment
        create_environment
    fi

    # Activate environment for remaining installations
    eval "$(conda shell.bash hook)"
    conda activate qiime2-amplicon-2024.2

    # Install dokdo from source
    install_dokdo

    # Verify imports
    verify_qiime2_imports

    # Test pipeline
    test_pipeline

    # Create scripts
    create_activation_script
    create_verification_script
    create_troubleshooting_script

    echo
    print_step "=================================================="
    print_status "🎉 INSTALLATION COMPLETED SUCCESSFULLY!"
    print_step "=================================================="
    echo
}

# Run main function
main