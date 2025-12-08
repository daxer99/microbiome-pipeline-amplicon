#!/usr/bin/env bash

set -e

echo "============================================"
echo " INSTALANDO QIIME2 2024.10 + DEBLUR FUNCIONAL"
echo "============================================"

# Detectar conda o mamba
if command -v mamba &> /dev/null; then
    MAMBA=mamba
    echo ">> Usando mamba"
else
    MAMBA=conda
    echo ">> Usando conda"
fi

ENV_NAME="qiime2-2024.10"

echo ">> Creando entorno $ENV_NAME"
$MAMBA create -y -n $ENV_NAME python=3.10

echo ">> Activando entorno"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

echo ">> Instalando QIIME2 2024.10"
$MAMBA install -y -c qiime2 -c conda-forge \
    qiime2=2024.10 \
    qiime-q2cli=2024.10 \
    q2templates=2024.10

echo ">> Instalando dependencias para Deblur"
$MAMBA install -y -c bioconda -c conda-forge \
    sortmerna=2.0 \
    vsearch \
    numpy \
    cython

echo ">> Instalando Deblur 1.1.1 desde bioconda"
$MAMBA install -y -c bioconda deblur=1.1.1

echo ">> Registrando plugin Deblur en QIIME2"
python - <<EOF
import qiime2.sdk
import q2_deblur
print("Deblur cargado correctamente:", q2_deblur)
EOF

echo "============================================"
echo "   QIIME2 + DEBLUR instalado correctamente"
echo "   Activar con: conda activate $ENV_NAME"
echo "============================================"
