#!/usr/bin/env bash

set -e

echo "==============================================="
echo " Instalando QIIME2 2024.2 + DEBLUR (FUNCIONAL) "
echo "==============================================="

ENV_NAME="qiime2-amplicon-2024.2"
YAML_URL="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2024.2/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"

echo ">> Eliminando entorno previo si existe..."
conda env remove -n $ENV_NAME -y > /dev/null 2>&1 || true

echo ">> Descargando YAML oficial..."
curl -s -L $YAML_URL -o qiime2.yml

echo ">> Creando entorno QIIME2 $ENV_NAME..."
mamba env create -n $ENV_NAME -f qiime2.yml

echo ">> Activando entorno..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

echo ">> Instalando Deblur desde Bioconda..."
mamba install -n $ENV_NAME -c bioconda -c conda-forge deblur -y

echo ">> Verificando import de Deblur..."
python3 - << 'EOF'
try:
    from qiime2.plugins.deblur.methods import denoise_16S
    print(">>> Deblur IMPORTADO correctamente dentro de QIIME2 ✔")
except Exception as e:
    print(">>> ERROR al importar Deblur ❌")
    print(e)
EOF

echo "==============================================="
echo " INSTALACIÓN COMPLETA 🎉"
echo " Para activar el entorno: "
echo "   conda activate $ENV_NAME"
echo "==============================================="
