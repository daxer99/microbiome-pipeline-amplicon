# Instalación manual paso a paso
echo "=== MANUAL INSTALLATION ==="

# 1. Create environment
echo "1. Creating environment..."
conda create -n qiime2-amplicon-2024.2 python=3.8 -y

# 2. Activate
echo "2. Activating..."
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2

# 3. Install QIIME2
echo "3. Installing QIIME2..."
conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released \
    qiime2 q2cli q2-demux q2-feature-table q2-taxa q2-deblur -y

# 4. Install PICRUSt2
echo "4. Installing PICRUSt2..."
conda install -c bioconda picrust2 vsearch cutadapt biom-format -y

# 5. Install Python packages
echo "5. Installing Python packages..."
conda install -c conda-forge pandas numpy matplotlib click biopython wget curl -y

# 6. Install dokdo
echo "6. Installing dokdo..."
cd /tmp
git clone https://github.com/sbslee/dokdo.git
cd dokdo
pip install .
cd ~/microbiome-pipeline-amplicon

# 7. Test
echo "7. Testing..."
python microbiome_cli.py --help

echo "=== DONE ==="