#!/usr/bin/env bash

set -e

echo "===================================="
echo "   CONFIGURANDO ENTORNO QIIME2"
echo "            Version 2025.10"
echo "===================================="

ENV_NAME="qiime2-2025.10"

echo ">> Detectando base de conda..."
source "$(conda info --base)/etc/profile.d/conda.sh"

echo ">> Creando entorno '${ENV_NAME}' con Python 3.11..."
mamba create -y -n ${ENV_NAME} python=3.11

echo ">> Activando entorno..."
conda activate ${ENV_NAME}

echo ">> Instalando QIIME2 2025.10 siguiendo el Quickstart oficial..."
mamba install -y -c qiime2 -c conda-forge \
    qiime2 q2cli q2templates

echo ">> Instalando plugins estándar del pipeline..."
mamba install -y -c qiime2 -c conda-forge \
    q2-dada2 \
    q2-vsearch \
    q2-diversity \
    q2-feature-classifier \
    q2-alignment \
    q2-phylogeny \
    q2-feature-table \
    q2-diversity \
    q2-metadata \
    q2-taxa \
    q2-types \
    q2-stats \
    q2-demux \
    q2-quality-control



echo ">> Instalando Cutadapt (requerido para trimming)..."
mamba install -y -c bioconda cutadapt

echo ">> Instalando dependencias externas útiles..."
mamba install -y -c conda-forge \
    pandas \
    biom-format \
    scikit-bio \
    matplotlib \
    wget

echo "===================================="
echo "     INSTALANDO DEBLUR MANUALMENTE"
echo "===================================="

echo ">> Instalando Deblur desde pip (método recomendado por los autores)..."
pip install deblur

echo ">> Instalando plugin de Deblur para QIIME2..."
pip install git+https://github.com/biocore/q2-deblur.git

echo ">> Registrando plugins..."
qiime dev refresh-cache

echo "===================================="
echo "  INSTALACIÓN COMPLETADA CORRECTAMENTE"
echo "===================================="
echo "Para usar QIIME2 2025.10 ejecuta:"
echo "       conda activate ${ENV_NAME}"
