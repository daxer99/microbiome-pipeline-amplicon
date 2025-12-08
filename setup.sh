#!/usr/bin/env bash
set -e

echo "===================================="
echo " INSTALANDO QIIME2 2025.10 + DEBLUR"
echo "    (compilado desde fuente)"
echo "===================================="

# Detect conda/mamba
if command -v mamba &> /dev/null; then
    CMD=mamba
else
    CMD=conda
fi

echo ">> Gestor detectado: $CMD"

ENV_NAME="qiime2-amplicon-2025.10"
YAML_URL="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"

echo ">> Creando entorno QIIME2 desde YAML oficial..."
$CMD env create -n $ENV_NAME --file $YAML_URL || echo ">> El entorno ya existe, continuando..."

echo ">> Activando entorno..."
eval "$($CMD shell hook)"
$CMD activate $ENV_NAME

echo "===================================="
echo " INSTALANDO DEPENDENCIAS PARA DEBLUR"
echo "===================================="

$CMD install -y -c conda-forge cython numpy h5py scipy

echo "===================================="
echo " COMPILANDO SORTMERNA 2.0"
echo "===================================="

cd /tmp
git clone https://github.com/biocore/sortmerna.git --branch v2.0 --depth 1
cd sortmerna
mkdir build && cd build
cmake ..
make -j4

echo ">> Instalando SortMeRNA en el entorno..."
cp sortmerna $CONDA_PREFIX/bin/
cp scripts/* $CONDA_PREFIX/bin/ 2>/dev/null || true

echo "===================================="
echo " COMPILANDO DEBLUR"
echo "===================================="

cd /tmp
git clone https://github.com/biocore/deblur.git --branch 1.1.1 --depth 1
cd deblur

$CONDA_PREFIX/bin/python -m pip install .

echo "===================================="
echo " VERIFICANDO INSTALACIÓN"
echo "===================================="

deblur --help && echo ">> Deblur instalado correctamente dentro del entorno" || echo ">> ERROR: Deblur no se instaló"

echo "===================================="
echo "   QIIME2 2025.10 + DEBLUR LISTO"
echo "===================================="
