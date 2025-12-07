
# Microbiome Pipeline Amplicon

A comprehensive pipeline for 16S rRNA microbiome analysis from raw sequencing data.

## Features

- **SRA Download**: Automated download of sequencing data from SRA.
- **Quality Control**: Comprehensive QC with FastQC and MultiQC.
- **Denoising**: ASV generation with Deblur.
- **Taxonomic Assignment**: Multiple classification methods.
- **Phylogenetic Analysis**: Tree building with MAFFT and FastTree.
- **Diversity Analysis**: Alpha and beta diversity metrics.
- **Metabolic Prediction**: PICRUSt2 for pathway inference.
- **Visualization**: Interactive tables, plots and reports.

### Prerequisites

- **Miniconda** or **Anaconda** installed.
- At least **10GB** of free disk space.
- O.S: Linux, macOS or Windows (WSL2 recommended).

## Installation

#### 1. Clone de repository

```bash
git clone https://github.com/daxer99/microbiome-pipeline-amplicon
cd microbiome-pipeline-amplicon
```
#### 2. Run the installation

```bash
chmod +x setup.sh
./setup.sh
```

#### 2.1. Manual installation (alternative)

```bash
# Create environment from YAML
conda env create -f environment.yml

# Activate environment
conda activate qiime2-amplicon-2024.2

# Install additional Python packages
pip install dokdo==1.16.0
```

#### 3. Install verification

```bash
# Activate environment
conda activate qiime2-amplicon-2024.2

# Test QIIME2
qiime --help

# Test PICRUSt2
picrust2_pipeline.py --version

# Test the pipeline
python microbiome_cli.py --help
```
## Data Resources

```
# Qiime2 documentation
https://docs.qiime2.org/2024.2/

# GreenGenes reference database
https://docs.qiime2.org/2024.2/data-resources/#greengenes-16s-rrna

# SILVA reference database
https://docs.qiime2.org/2024.2/data-resources/#silva-16s-18s-rrna

# UNITE reference database
https://docs.qiime2.org/2024.2/data-resources/#unite-fungal-its

# EzBioCloud reference database
https://www.ezbiocloud.net/resources/16s_download
```
## Usage

```python
# 1. Download sequences from SRA
python microbiome_cli.py download example/sra.csv --output-dir example/samples

# 2. Create manifest file
python microbiome_cli.py create-manifest --input-dir example/samples --output-file example/manifest.csv

# 3. Import to QIIME2
python microbiome_cli.py import-sample-seqs example/manifest.csv --output-dir example/demux.qza

# 4. Quality control
python microbiome_cli.py quality-control example/demux.qza --output-dir example/qc

# 5. Denoising with Deblur
python microbiome_cli.py run-deblur example/demux.qza --left-trim-len 10 --trim-length 125 --jobs-to-start 8 --output-dir example/deblur_results

# 5.1. Alternative
python microbiome_cli.py run-deblur example/qc/flitered_seqs.qza --left-trim-len 0 --trim-length 250 --jobs-to-start 8 --output-dir example/deblur_results

# 6. Phylogenetic tree
python microbiome_cli.py build-phylogeny example/deblur_results/rep-seqs.qza --output-dir example/phylogeny

# 7. Taxonomic assignment
python microbiome_cli.py assign-taxonomy example/deblur_results/table.qza example/deblur_results/rep-seqs.qza ref_db/ref-seqs.qza ref_db/ref-taxa.qza example/metadata.tsv --cpus 8 --output-dir example/taxa

# 8.1 Alpha Diversity analysis
python microbiome_cli.py alpha-diversity example/deblur_results/table.qza --metrics ace,gini_index,shannon --output-dir example/alpha

python microbiome_cli.py alpha-diversity example/deblur_results/table.qza --metrics ace,faith_pd,shannon --rooted-tree example/phylogeny/rooted_tree.qza --output-dir example/alpha

# 8.2 Beta Diversity analysis
python microbiome_cli.py beta-diversity example/deblur_results/table.qza --metrics braycurtis --metadata example/metadata.tsv --hue subject --output-dir example/beta

python microbiome_cli.py beta-diversity example/deblur_results/table.qza --metrics braycurtis  --phylo-metrics unweighted_unifrac --rooted-tree example/phylogeny/rooted_tree.qza --metadata example/metadata.tsv --hue subject --output-dir example/beta

# 9. Metabolic pathway prediction
python microbiome_cli.py predict-metabolic-pathways example/deblur_results/table.qza example/deblur_results/rep-seqs.qza --threads 8 --output-dir example/picrust2
```


## Citation


If you use this pipeline, please cite:

- **QIIME2:** Bolyen, E. et al. (2019) “Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2,” Nature Biotechnology, 37(8), pp. 852–857.

- **Deblur:** Amir, A. et al. (2017) “Deblur Rapidly Resolves Single-Nucleotide Community Sequence Patterns,” mSystems, 2(2), p. 10.1128/msystems.00191-16.

- **PICRUSt2:** Douglas, G.M., et al. PICRUSt2 for prediction of metagenome functions. Nat Biotechnol 38, 685–688 (2020).

## License

[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)


## Authors

- [@Rodrigo Peralta](https://www.linkedin.com/in/rodrigo-peralta-28401212b/)

