#!/usr/bin/env bash
set -e

echo "===================================="
echo " CONFIGURANDO ENTORNO QIIME2 2025.10"
echo "===================================="

ENV_NAME="qiime2-2025.10"

# Detectar gestor: micromamba → mamba → conda
if command -v micromamba &>/dev/null; then
    MANAGER="micromamba"
    RUN="micromamba run -n $ENV_NAME"
elif command -v mamba &>/dev/null; then
    MANAGER="mamba"
    RUN="mamba run -n $ENV_NAME"
else
    MANAGER="conda"
    RUN="conda run -n $ENV_NAME"
fi

echo ">> Gestor detectado: $MANAGER"
echo ">> Creando entorno $ENV_NAME..."

$MANAGER create -y -n $ENV_NAME python=3.11

echo ">> Instalando QIIME2 2025.10..."
$MANAGER install -y -n $ENV_NAME -c qiime2 -c conda-forge \
    qiime2 q2cli q2templates

echo ">> Instalando plugins comunes..."
$MANAGER install -y -n $ENV_NAME -c qiime2 -c conda-forge \
    q2-diversity q2-metadata q2-feature-classifier \
    q2-alignment q2-phylogeny q2-vsearch q2-dada2

echo ">> Instalando Cutadapt..."
$MANAGER install -y -n $ENV_NAME -c bioconda cutadapt

echo "===================================="
echo " INSTALANDO DEBLUR"
echo "===================================="

# pip dentro del entorno SIN activarlo
$RUN pip install deblur
$RUN pip install git+https://github.com/biocore/q2-deblur.git

# reconstruir cache de QIIME2
$RUN qiime dev refresh-cache

echo "===================================="
echo " INSTALACIÓN COMPLETA"
echo ""
echo " Para activar el entorno:"
echo ""
echo "   $MANAGER activate $ENV_NAME"
echo ""
echo "===================================="
