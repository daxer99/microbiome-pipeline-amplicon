#!/usr/bin/env bash
set -e

echo "==============================="
echo " Instalando QIIME2 2024.10"
echo "==============================="

# Cargar conda
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

ENV="qiime2-2024.10"

echo ">> Eliminando entorno previo si existe..."
conda env remove -n $ENV -y || true

echo ">> Creando entorno desde el YAML oficial..."
conda env create \
  -n $ENV \
  -f https://raw.githubusercontent.com/qiime2/distributions/2024.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml

echo ">> Activando entorno..."
conda activate $ENV

echo ">> Verificando instalación de qiime CLI..."
qiime --help >/dev/null 2>&1 && echo "QIIME CLI OK" || { echo "ERROR: qiime no está instalado"; exit 1; }

echo ">> Instalando Deblur dentro del entorno..."
conda install -y -c bioconda -c conda-forge deblur=1.1.1

echo ">> Verificando plugin Deblur..."
python - <<EOF
from qiime2.plugins.deblur.methods import denoise_16S
print("Deblur OK")
EOF

echo "=================================="
echo " QIIME2 + Deblur instalados OK"
echo " Activar con: conda activate $ENV"
echo "=================================="
