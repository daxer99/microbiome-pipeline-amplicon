#!/bin/bash

# ==============================================================================
# Setup script for Microbiome Pipeline Amplicon - QIIME2 2024.2
# Official installation method + channel_priority flexible fix
# Compatible with Linux and WSL2
# ==============================================================================

set -e  # Exit on error

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Configuration
ENV_NAME="qiime2-amplicon-2024.2"
QIIME_VERSION="2024.2"
QIIME_YML_URL="https://raw.githubusercontent.com/qiime2/distributions/dev/${QIIME_VERSION}/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"

clear
echo -e "${CYAN}=============================================="
echo "   Microbiome Pipeline Amplicon"
echo "   QIIME2 ${QIIME_VERSION} Installation"
echo "=============================================="
echo -e "${NC}"

# Check if running on Linux/WSL
if [[ ! "$(uname -s)" == "Linux" ]]; then
    echo -e "${YELLOW}⚠ Warning: This script is optimized for Linux/WSL2${NC}"
    echo "For macOS, use the official documentation:"
    echo "https://library.qiime2.org/quickstart/amplicon"
    echo ""
    read -p "Continue anyway? (y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 0
    fi
fi

# Check conda
if ! command -v conda &> /dev/null; then
    echo -e "${RED}❌ Error: conda not found${NC}"
    echo ""
    echo "Please install Miniconda first:"
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh"
    echo ""
    exit 1
fi

echo -e "${GREEN}✓${NC} conda found: $(conda --version)"

# Update conda
echo ""
echo "Updating conda..."
conda update -n base -c defaults conda -y 2>&1 | tail -3

# CRITICAL: Set channel priority to flexible
echo ""
echo -e "${CYAN}=============================================="
echo "Configuring conda channels..."
echo "=============================================="
echo -e "${NC}"

echo -e "${BLUE}Setting channel_priority to flexible...${NC}"
echo "This is REQUIRED to resolve QIIME2 2024.2 dependencies"
conda config --set channel_priority flexible

# Verify configuration
PRIORITY=$(conda config --show channel_priority | grep -o 'flexible\|strict')
if [[ "$PRIORITY" == "flexible" ]]; then
    echo -e "${GREEN}✓ channel_priority: flexible${NC}"
else
    echo -e "${RED}❌ Failed to set channel_priority${NC}"
    exit 1
fi

# Remove old environment if exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo ""
    echo -e "${YELLOW}⚠ Environment '${ENV_NAME}' already exists${NC}"
    read -p "Remove and recreate? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing old environment..."
        conda deactivate 2>/dev/null || true
        conda env remove -n ${ENV_NAME} -y
    else
        echo "Installation cancelled."
        exit 0
    fi
fi

# Clean conda cache
echo ""
echo "Cleaning conda cache..."
conda clean --all -y > /dev/null 2>&1

# Create environment using official method
echo ""
echo -e "${CYAN}=============================================="
echo "Creating QIIME2 environment..."
echo "Using official installation method"
echo "This may take 15-25 minutes..."
echo "=============================================="
echo -e "${NC}"
echo ""
echo "Command:"
echo -e "${BLUE}conda env create \\"
echo "  --name ${ENV_NAME} \\"
echo "  --file ${QIIME_YML_URL}${NC}"
echo ""

if conda env create \
    --name ${ENV_NAME} \
    --file ${QIIME_YML_URL}; then
    echo ""
    echo -e "${GREEN}✓ QIIME2 environment created successfully${NC}"
else
    echo ""
    echo -e "${RED}❌ Failed to create QIIME2 environment${NC}"
    echo ""
    echo "Troubleshooting:"
    echo "1. Check internet connection"
    echo "2. Verify channel_priority: conda config --show channel_priority"
    echo "3. Try: conda clean --all && ./setup.sh"
    echo "4. Check QIIME2 forum: https://forum.qiime2.org/"
    exit 1
fi

# Activate environment
echo ""
echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${ENV_NAME}

if [[ $CONDA_DEFAULT_ENV != ${ENV_NAME} ]]; then
    echo -e "${RED}❌ Failed to activate environment${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Environment activated${NC}"

# Install additional tools
echo ""
echo -e "${CYAN}=============================================="
echo "Installing additional tools..."
echo "=============================================="
echo -e "${NC}"

# PICRUSt2
echo ""
echo "Installing PICRUSt2..."
if conda install -c bioconda -c conda-forge picrust2 -y 2>&1 | tail -5; then
    if command -v picrust2_pipeline.py &> /dev/null; then
        echo -e "${GREEN}✓ PICRUSt2 installed${NC}"
    else
        echo -e "${YELLOW}⚠ PICRUSt2 installation completed but not in PATH${NC}"
    fi
else
    echo -e "${YELLOW}⚠ PICRUSt2 installation had issues (optional)${NC}"
fi

# SRA Toolkit
echo ""
echo "Installing SRA Toolkit..."
if conda install -c bioconda sra-tools -y 2>&1 | tail -5; then
    if command -v fastq-dump &> /dev/null; then
        echo -e "${GREEN}✓ SRA Toolkit installed${NC}"
    else
        echo -e "${YELLOW}⚠ SRA Toolkit installation completed but not in PATH${NC}"
    fi
else
    echo -e "${YELLOW}⚠ SRA Toolkit installation had issues (optional)${NC}"
fi

# Python packages
echo ""
echo "Installing Python packages..."
echo -n "  dokdo... "
if pip install --quiet dokdo==1.16.0 2>/dev/null; then
    echo -e "${GREEN}✓${NC}"
else
    echo -e "${YELLOW}⚠${NC}"
fi

# Verification
echo ""
echo -e "${CYAN}=============================================="
echo "Verifying installation..."
echo "=============================================="
echo -e "${NC}"

# Test QIIME2
echo ""
if qiime --version &> /dev/null; then
    QIIME_V=$(qiime --version 2>&1)
    echo -e "${GREEN}✓ QIIME2: ${QIIME_V}${NC}"
else
    echo -e "${RED}❌ QIIME2 verification failed${NC}"
    exit 1
fi

# Test plugins
echo ""
echo "Checking QIIME2 plugins:"
PLUGINS_OK=0
PLUGINS_TOTAL=0

declare -a plugins=("deblur" "dada2" "feature-table" "diversity" "taxa" "phylogeny" "demux")

for plugin in "${plugins[@]}"; do
    PLUGINS_TOTAL=$((PLUGINS_TOTAL + 1))
    if qiime ${plugin} --help &> /dev/null; then
        echo -e "  ${GREEN}✓${NC} q2-${plugin}"
        PLUGINS_OK=$((PLUGINS_OK + 1))
    else
        echo -e "  ${RED}✗${NC} q2-${plugin}"
    fi
done

echo ""
echo -e "Plugins working: ${PLUGINS_OK}/${PLUGINS_TOTAL}"

# Test optional tools
echo ""
echo "Optional tools:"
if command -v picrust2_pipeline.py &> /dev/null; then
    PICRUST_V=$(picrust2_pipeline.py --version 2>&1 | head -n 1 | cut -d' ' -f2)
    echo -e "  ${GREEN}✓${NC} PICRUSt2 ${PICRUST_V}"
else
    echo -e "  ${YELLOW}⚠${NC} PICRUSt2 not available"
fi

if command -v fastq-dump &> /dev/null; then
    SRA_V=$(fastq-dump --version 2>&1 | head -n 1 | cut -d' ' -f3)
    echo -e "  ${GREEN}✓${NC} SRA Toolkit ${SRA_V}"
else
    echo -e "  ${YELLOW}⚠${NC} SRA Toolkit not available"
fi

# Test microbiome_cli.py
echo ""
if [[ -f "microbiome_cli.py" ]]; then
    if python microbiome_cli.py --help &> /dev/null; then
        echo -e "${GREEN}✓${NC} microbiome_cli.py is functional"
    else
        echo -e "${YELLOW}⚠${NC} microbiome_cli.py has issues"
    fi
else
    echo -e "${YELLOW}⚠${NC} microbiome_cli.py not found in current directory"
fi

# Success summary
echo ""
echo -e "${CYAN}=============================================="
echo -e "${GREEN}✓ Installation completed successfully!${NC}"
echo -e "${CYAN}=============================================="
echo -e "${NC}"

# Create activation script
cat > activate_qiime2.sh << 'ACTIVATE_SCRIPT'
#!/bin/bash
# Quick activation script for QIIME2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2024.2
echo "✓ QIIME2 environment activated"
echo "To test: qiime --version"
ACTIVATE_SCRIPT

chmod +x activate_qiime2.sh
echo -e "${GREEN}✓${NC} Created activation script: ${BLUE}./activate_qiime2.sh${NC}"

# Usage instructions
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${CYAN}Quick Start Guide${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo -e "${BLUE}1. Activate environment:${NC}"
echo -e "   ${GREEN}conda activate ${ENV_NAME}${NC}"
echo -e "   or: ${GREEN}./activate_qiime2.sh${NC}"
echo ""
echo -e "${BLUE}2. Test installation:${NC}"
echo -e "   ${GREEN}qiime info${NC}"
echo ""
echo -e "${BLUE}3. Run your pipeline:${NC}"
echo -e "   ${GREEN}python microbiome_cli.py --help${NC}"
echo ""
echo -e "${BLUE}4. Deactivate when done:${NC}"
echo -e "   ${GREEN}conda deactivate${NC}"
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${CYAN}Environment Details${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo -e "  Name: ${GREEN}${ENV_NAME}${NC}"
echo -e "  Python: ${GREEN}$(python --version 2>&1)${NC}"
echo -e "  Location: ${GREEN}$(conda info --base)/envs/${ENV_NAME}${NC}"
echo -e "  Channel priority: ${GREEN}flexible${NC} (required for QIIME2 2024.2)"
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${CYAN}Resources${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo "  • QIIME2 Documentation:"
echo "    https://library.qiime2.org/quickstart/amplicon"
echo ""
echo "  • Your Repository:"
echo "    https://github.com/daxer99/microbiome-pipeline-amplicon"
echo ""
echo "  • QIIME2 Forum (for help):"
echo "    https://forum.qiime2.org/"
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""