#!/usr/bin/env bash
set -e

echo "============================================"
echo " INSTALANDO QIIME2 2024.10 + DEBLUR FUNCIONAL"
echo "============================================"

# Forzar siempre conda (micromamba rompe el entorno)
CONDA=$(conda info --base)
source $CONDA/etc/profile.d/conda.sh

MAMBA="conda"   # <--- fuerza el uso de conda SIEMPRE

ENV="qiime2-2024.10"

echo ">> Creando entorno con CONDA"
conda create -y -n $ENV python=3.10

echo ">> Activando entorno"
conda activate $ENV

echo ">> Instalando QIIME2 2024.10"
conda install -y -c qiime2 -c conda-forge \
    qiime2=2024.10 \
    qiime-q2cli=2024.10 \
    q2templates=2024.10

echo ">> Instalando dependencias Deblur"
conda install -y -c bioconda -c conda-forge \
    sortmerna=2.0 \
    vsearch \
    numpy \
    cython

echo ">> Instalando Deblur"
conda install -y -c bioconda deblur=1.1.1

echo ">> Verificando el plugin"
python - <<EOF
from qiime2.plugins.deblur.methods import denoise_16S
print("Deblur cargado OK")
EOF

echo "============================================"
echo "     QIIME2 + Deblur instalados exitosamente"
echo "     Ejecutar: conda activate $ENV"
echo "============================================"
