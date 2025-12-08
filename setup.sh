#!/usr/bin/env bash
set -e

echo "==============================="
echo " Instalando QIIME2 2024.10"
echo "==============================="

# Cargar conda
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

ENV="qiime2-2024.10"
YAML_URL="https://raw.githubusercontent.com/qiime2/distributions/2024.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"

echo ">> Eliminando entorno previo si existe..."
conda env remove -n $ENV -y || true

echo ">> Descargando YAML oficial..."
curl -fsSL "$YAML_URL" -o qiime2.yml

echo ">> Creando entorno desde YAML..."
# NOTA: AHORA ES OBLIGATORIO AGREGAR:
#   --environment-spec environment.yml
mamba env create \
    --name $ENV \
    --file qiime2.yml \
    --environment-spec environment.yml

echo ">> Activando entorno..."
conda activate $ENV

echo ">> Verificando 'qiime'..."
qiime --help >/dev/null 2>&1 && echo "QIIME CLI OK" || { echo "ERROR: qiime no está instalado"; exit 1; }

echo ">> Instalando Deblur..."
mamba install -y -c bioconda -c conda-forge deblur=1.1.1

echo ">> Probando import de Deblur..."
python - <<EOF
from qiime2.plugins.deblur.methods import denoise_16S
print("Deblur OK ✔")
EOF

echo "====================================="
echo " QIIME2 + Deblur instalados correctamente."
echo " Activar con: conda activate $ENV"
echo "====================================="
