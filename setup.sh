#!/usr/bin/env bash

set -e

echo "======================================"
echo "  Creando entorno QIIME2 2024.2"
echo "======================================"

# Crear archivo YAML
cat << 'EOF' > qiime2-2024.2.yml
name: qiime2-2024.2
channels:
  - qiime2/label/r2024.2
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - qiime2-core=2024.2
  - q2cli=2024.2
  - q2templates=2024.2
  - qiime2=2024.2

  - q2-dada2=2024.2
  - q2-demux=2024.2
  - q2-feature-classifier=2024.2
  - q2-metadata=2024.2
  - q2-diversity=2024.2
  - q2-emperor=2024.2
  - q2-taxa=2024.2
  - q2-types=2024.2

  - python=3.10
  - numpy
  - scipy
  - pandas
  - scikit-bio
  - biom-format
  - statsmodels
  - matplotlib
  - seaborn
  - jupyterlab
  - ipython

  - fasttree
  - mafft
  - vsearch

  - wget
  - pip

  - pip:
    - biopython
    - plotnine
EOF

echo "[OK] Archivo qiime2-2024.2.yml generado."

echo "======================================"
echo "  Creando entorno..."
echo "======================================"

conda env create -f qiime2-2024.2.yml

echo "======================================"
echo "  Activando entorno..."
echo "======================================"

# Activación según shell
if [[ $SHELL == *"bash"* ]]; then
    echo "Detectado bash."
    source ~/.bashrc || true
elif [[ $SHELL == *"zsh"* ]]; then
    echo "Detectado zsh."
    source ~/.zshrc || true
fi

conda activate qiime2-2024.2

echo "======================================"
echo "  Verificando instalación"
echo "======================================"

qiime --help || { echo "ERROR: QIIME2 no se instaló correctamente"; exit 1; }

echo "======================================"
echo "  QIIME2 2024.2 instalado correctamente"
echo "======================================"
echo "Entorno activo: qiime2-2024.2"
