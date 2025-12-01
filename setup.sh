#!/bin/bash

# Microbiome Pipeline Amplicon - Setup Script
# Optimized for WSL2 and automatic error recovery

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

# Create Conda environment with retry logic
create_environment() {
    local env_name="qiime2-amplicon-2024.2"
    local max_retries=3
    local retry_count=0

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

    # Try to create environment with retries
    while [ $retry_count -lt $max_retries ]; do
        print_status "Attempt $((retry_count + 1)) of $max_retries..."

        if [ -f "environment.yml" ]; then
            # First, try with mamba if available (much faster)
            if command -v mamba &> /dev/null; then
                print_status "Using mamba (faster package manager)..."
                if mamba env create -f environment.yml; then
                    print_status "Environment created successfully with mamba"
                    return 0
                fi
            else
                # Try to install mamba
                print_status "Installing mamba for faster package management..."
                conda install -c conda-forge mamba -y
                if mamba env create -f environment.yml; then
                    print_status "Environment created successfully with mamba"
                    return 0
                fi
            fi

            # If mamba fails, try with conda
            print_status "Falling back to conda..."
            if conda env create -f environment.yml; then
                print_status "Environment created successfully with conda"
                return 0
            fi
        else
            print_error "environment.yml not found in current directory"
            exit 1
        fi

        retry_count=$((retry_count + 1))
        if [ $retry_count -lt $max_retries ]; then
            print_warning "Creation failed. Retrying in 10 seconds..."
            sleep 10
        fi
    done

    print_error "Failed to create environment after $max_retries attempts"

    # Try alternative simplified environment
    print_status "Trying alternative simplified installation..."
    create_simplified_environment
}

# Create simplified environment if full one fails
create_simplified_environment() {
    local env_name="qiime2-amplicon-2024.2"

    print_status "Creating simplified environment..."

    # Create basic environment
    conda create -n "$env_name" python=3.8 -y

    # Activate environment
    eval "$(conda shell.bash hook)"
    conda activate "$env_name"

    # Install core packages individually
    print_status "Installing core packages..."

    # QIIME2 core
    conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released qiime2 q2cli q2-feature-table q2-taxa -y

    # PICRUSt2
    conda install -c bioconda picrust2 -y

    # Essential bioinformatics tools
    conda install -c bioconda vsearch cutadapt fastqc -y

    # Python packages
    conda install -c conda-forge pandas numpy matplotlib seaborn click -y

    print_status "Simplified environment created successfully"
    print_warning "Note: Some advanced features may not be available"
}

# Install additional Python packages
install_python_packages() {
    print_status "Installing additional Python packages..."

    # Activate environment
    eval "$(conda shell.bash hook)"
    if conda activate qiime2-amplicon-2024.2; then
        pip install dokdo==1.16.0
        print_status "Python packages installed successfully"
    else
        print_warning "Could not activate environment for pip installation"
        print_status "You can install dokdo later with: pip install dokdo==1.16.0"
    fi
}

# Verify installation
verify_installation() {
    print_status "Verifying installation..."

    # Activate environment first
    eval "$(conda shell.bash hook)"
    if ! conda activate qiime2-amplicon-2024.2; then
        print_error "Failed to activate environment for verification"
        return 1
    fi

    local success=true

    # Check QIIME2
    if qiime --help &> /dev/null; then
        print_status "✓ QIIME2 installed correctly"
    else
        print_error "✗ QIIME2 installation failed"
        success=false
    fi

    # Check PICRUSt2
    if picrust2_pipeline.py --version &> /dev/null; then
        print_status "✓ PICRUSt2 installed correctly"
    else
        print_error "✗ PICRUSt2 installation failed"
        success=false
    fi

    # Check other key tools
    tools=("fastqc" "vsearch" "cutadapt")
    for tool in "${tools[@]}"; do
        if $tool --version &> /dev/null; then
            print_status "✓ $tool installed correctly"
        else
            print_warning "⚠ $tool version check failed (but may still work)"
        fi
    done

    if $success; then
        print_status "✓ Core verification passed"
    else
        print_warning "⚠ Some verifications failed, but pipeline may still work"
    fi
}

# Create activation script
create_activation_script() {
    cat > activate_pipeline.sh << 'EOF'
#!/bin/bash
# Script to activate the microbiome pipeline environment

echo "Activating Microbiome Pipeline Environment..."

# Initialize conda
eval "$(conda shell.bash hook)"

# Try to activate environment
if conda activate qiime2-amplicon-2024.2; then
    echo "=================================================="
    echo "Environment activated successfully!"
    echo "=================================================="
    echo ""
    echo "Available commands:"
    echo "  python microbiome_cli.py --help"
    echo "  qiime --help"
    echo "  picrust2_pipeline.py --version"
    echo ""
    echo "Example workflow:"
    echo "  1. python microbiome_cli.py download samples.csv"
    echo "  2. python microbiome_cli.py quality-control demux.qza"
    echo "  3. python microbiome_cli.py predict-metabolic-pathways table.qza rep-seqs.qza"
else
    echo "Failed to activate environment."
    echo "Try: conda activate qiime2-amplicon-2024.2"
fi
EOF

    chmod +x activate_pipeline.sh
    print_status "Created activation script: ./activate_pipeline.sh"
}

# Create troubleshooting script
create_troubleshooting_script() {
    cat > troubleshoot.sh << 'EOF'
#!/bin/bash
echo "Microbiome Pipeline Troubleshooting"
echo "==================================="

# Check conda
echo "1. Checking Conda..."
if command -v conda &> /dev/null; then
    echo "   ✓ Conda is installed"
else
    echo "   ✗ Conda not found"
fi

# Check environment
echo "2. Checking environment..."
if conda env list | grep -q "qiime2-amplicon-2024.2"; then
    echo "   ✓ Environment exists"

    # Activate and check tools
    eval "$(conda shell.bash hook)"
    if conda activate qiime2-amplicon-2024.2; then
        echo "3. Checking tools in environment:"

        # List of tools to check
        declare -A tools
        tools=(
            ["QIIME2"]="qiime --help"
            ["PICRUSt2"]="picrust2_pipeline.py --version"
            ["vsearch"]="vsearch --version"
            ["cutadapt"]="cutadapt --version"
            ["Python CLI"]="python microbiome_cli.py --help"
        )

        for tool_name in "${!tools[@]}"; do
            command_str="${tools[$tool_name]}"
            if eval "$command_str" &> /dev/null; then
                echo "   ✓ $tool_name works"
            else
                echo "   ✗ $tool_name NOT working"
            fi
        done
    else
        echo "   ✗ Could not activate environment"
    fi
else
    echo "   ✗ Environment not found"
fi

echo ""
echo "Common solutions:"
echo "1. If environment missing: ./setup.sh"
echo "2. If tools missing: conda install -c bioconda <tool>"
echo "3. For disk space: conda clean --all"
EOF

    chmod +x troubleshoot.sh
    print_status "Created troubleshooting script: ./troubleshoot.sh"
}

# Main installation function
main() {
    echo
    print_status "Starting installation process..."

    # Check conda
    check_conda

    # Clean conda cache to free space
    print_status "Cleaning conda cache..."
    conda clean --all -y

    # Create environment
    create_environment

    # Install additional packages
    install_python_packages

    # Verify installation
    verify_installation

    # Create activation script
    create_activation_script

    # Create troubleshooting script
    create_troubleshooting_script

    echo
    print_status "=================================================="
    print_status "Installation completed! 🎉"
    print_status "=================================================="
    echo
    print_status "NEXT STEPS:"
    print_status "1. Always start with: source ./activate_pipeline.sh"
    print_status "2. Test installation: ./troubleshoot.sh"
    print_status "3. View commands: python microbiome_cli.py --help"
    echo
    print_status "For WSL2 users:"
    print_status "- Store data in /mnt/c/ for Windows access"
    print_status "- Use: \\\\wsl$\\Ubuntu\\home\\$(whoami) in Windows Explorer"
    echo
}

# Run main function
main