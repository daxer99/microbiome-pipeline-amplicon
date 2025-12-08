#!/usr/bin/env bash

set -e

ENV_NAME="qiime2-amplicon-2024.2"
YAML_URL="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2024.2/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"

echo "==============================================="
echo " Instalación de QIIME2 2024.2 (Amplicon)"
echo "==============================================="
echo ""

# Verificar conda
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda no está disponible en el PATH."
    exit 1
fi

# Usar mamba si existe
if command -v mamba &> /dev/null; then
    CMD=mamba
    echo ">> Usando mamba"
else
    CMD=conda
    echo ">> Usando conda"
fi

echo ""
echo ">> Ajustando channel_priority = flexible"
conda config --set channel_priority flexible

echo ""
echo ">> Eliminando entorno previo (si existe)..."
conda remove -n "$ENV_NAME" --all -y || true

echo ""
echo ">> Creando entorno desde el YAML oficial..."
$CMD env create \
    --name "$ENV_NAME" \
    --file "$YAML_URL"

echo ""
echo ">> Activando entorno para verificar instalación..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo ""
echo ">> Verificando instalación de QIIME2..."
if qiime --help &> /dev/null; then
    echo "==============================================="
    echo " QIIME2 2024.2 instalado correctamente 🎉"
    echo " Entorno: $ENV_NAME"
    echo " Para activarlo: conda activate $ENV_NAME"
    echo "==============================================="
else
    echo "ERROR: QIIME2 no se instaló correctamente."
    exit 1
fi
