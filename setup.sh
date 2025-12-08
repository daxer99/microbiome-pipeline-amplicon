#!/usr/bin/env bash
set -e

echo "===================================="
echo " CONFIGURANDO ENTORNO QIIME2 2025.10"
echo "===================================="

ENV_NAME="qiime2-2025.10"

# Detectar gestor (micromamba → mamba → conda)
if command -v micromamba &> /dev/null; then
    MANAGER="micromamba"
    ACTIVATE="micromamba activate"
elif command -v mamba &> /dev/null; then
    MANAGER="mamba"
    ACTIVATE="mamba activate"
else
    MANAGER="conda"
    ACTIVATE="conda activate"
fi

echo ">> Gestor detectado: $MANAGER"
echo ">> Creando entorno $ENV_NAME..."

$MANAGER create -y -n $ENV_NAME python=3.11

echo ">> Activando entorno..."
source "$($MANAGER info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
$ACTIVATE $ENV_NAME

echo ">> Instalando QIIME2 2025.10..."
$MANAGER install -y -c qiime2 -c conda-forge \
    qiime2 q2cli q2templates

echo ">> Instalando plugins comunes..."
$MANAGER install -y -c qiime2 -c conda-forge \
    q2-diversity q2-metadata q2-feature-classifier \
    q2-alignment q2-phylogeny q2-vsearch q2-dada2

echo ">> Instalando Cutadapt..."
$MANAGER install -y -c bioconda cutadapt

echo "===================================="
echo "       INSTALANDO DEBLUR"
echo "===================================="

pip install deblur
pip install git+https://github.com/biocore/q2-deblur.git

qiime dev refresh-cache

echo "===================================="
echo " INSTALACIÓN COMPLETA"
echo " Activa tu entorno con:"
echo "     $ACTIVATE $ENV_NAME"
echo "===================================="
