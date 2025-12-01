#!/bin/bash

# Microbiome Pipeline Amplicon - Setup Script
# This script sets up the complete environment for the 16S microbiome analysis pipeline

set -e  # Exit on any error

echo "=================================================="
echo "Microbiome Pipeline Amplicon - Installation Script"
echo "=================================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
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

# Check if Conda is installed
check_conda() {
    if ! command -v conda &> /dev/null; then
        print_error "Conda is not installed. Please install Miniconda or Anaconda first."
        echo "Download Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
    print_status "Conda is installed"
}

# Create Conda environment
create_environment() {
    local env_name="qiime2-amplicon-2024.2"

    print_status "Creating Conda environment: $env_name"

    # Check if environment already exists
    if conda env list | grep -q "$env_name"; then
        print_warning "Environment $env_name already exists"
        read -p "Do you want to remove and recreate it? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            conda env remove -n "$env_name"
        else
            print_status "Using existing environment"
            return 0
        fi
    fi

    # Create environment from YAML file
    if [ -f "environment.yml" ]; then
        conda env create -f environment.yml
    else
        print_error "environment.yml not found in current directory"
        exit 1
    fi

    print_status "Environment created successfully"
}

# Install additional Python packages
install_python_packages() {
    print_status "Installing additional Python packages..."

    # Activate environment and install packages
    conda run -n qiime2-amplicon-2024.2 pip install dokdo==1.16.0

    print_status "Python packages installed successfully"
}

# Verify installation
verify_installation() {
    print_status "Verifying installation..."

    # Check QIIME2
    if conda run -n qiime2-amplicon-2024.2 qiime --help &> /dev/null; then
        print_status "✓ QIIME2 installed correctly"
    else
        print_error "✗ QIIME2 installation failed"
        exit 1
    fi

    # Check PICRUSt2
    if conda run -n qiime2-amplicon-2024.2 picrust2_pipeline.py --version &> /dev/null; then
        print_status "✓ PICRUSt2 installed correctly"
    else
        print_error "✗ PICRUSt2 installation failed"
        exit 1
    fi

    # Check other key tools
    tools=("fastqc" "vsearch" "cutadapt" "mafft" "fasttree")
    for tool in "${tools[@]}"; do
        if conda run -n qiime2-amplicon-2024.2 $tool --version &> /dev/null; then
            print_status "✓ $tool installed correctly"
        else
            print_warning "⚠ $tool version check failed (but may still work)"
        fi
    done

    print_status "All verifications completed"
}

# Create activation script
create_activation_script() {
    cat > activate_pipeline.sh << 'EOF'
#!/bin/bash
# Script to activate the microbiome pipeline environment

echo "Activating Microbiome Pipeline Environment..."
conda activate qiime2-amplicon-2024.2

# Check if activation was successful
if [ $? -eq 0 ]; then
    echo "Environment activated successfully!"
    echo "You can now run: python microbiome_cli.py --help"
else
    echo "Failed to activate environment. Please check the installation."
fi
EOF

    chmod +x activate_pipeline.sh
    print_status "Created activation script: ./activate_pipeline.sh"
}

# Main installation function
main() {
    echo
    print_status "Starting installation process..."

    # Check conda
    check_conda

    # Create environment
    create_environment

    # Install additional packages
    install_python_packages

    # Verify installation
    verify_installation

    # Create activation script
    create_activation_script

    echo
    print_status "=================================================="
    print_status "Installation completed successfully! 🎉"
    print_status "=================================================="
    echo
    print_status "To use the pipeline:"
    print_status "1. Activate the environment: ./activate_pipeline.sh"
    print_status "2. Or manually: conda activate qiime2-amplicon-2024.2"
    print_status "3. Run: python microbiome_cli.py --help"
    echo
    print_status "Available commands:"
    print_status "  download, quality-control, run-deblur, assign-taxonomy"
    print_status "  build-phylogeny, alpha-diversity, beta-diversity"
    print_status "  predict-metabolic-pathways"
    echo
}

# Run main function
main