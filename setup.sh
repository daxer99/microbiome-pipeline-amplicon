#!/usr/bin/env bash
set -e

echo "===================================="
echo " CONFIGURANDO QIIME2 2025.10 + DEBLUR"
echo "===================================="

# Detectar gestor
if command -v mamba &> /dev/null; then
    PM="mamba"
    echo ">> Gestor detectado: mamba"
else
    PM="conda"
    echo ">> Gestor detectado: conda"
fi

ENV_NAME="qiime2-amplicon-2025.10"
YAML_URL="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"

# Crear entorno si no existe
if ! conda env list | grep -q "$ENV_NAME"; then
    echo ">> Creando entorno $ENV_NAME usando YAML oficial..."
    $PM env create --name $ENV_NAME --file "$YAML_URL"
else
    echo ">> El entorno $ENV_NAME ya existe. Continuando..."
fi

echo ">> Activando entorno..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

echo ">> Instalando Deblur..."
$PM install -y -c bioconda deblur

echo ">> Instalación completa!"
echo "------------------------------------"
echo " Para usar QIIME2 ejecutar:"
echo "   conda activate $ENV_NAME"
echo "------------------------------------"
