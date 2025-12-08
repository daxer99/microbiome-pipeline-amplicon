#!/usr/bin/env bash
set -e

echo "==========================================="
echo " Instalación de QIIME2 2024.5 + Deblur"
echo "==========================================="

ENV_NAME="qiime2-amplicon-2024.5"
YAML_URL="https://packages.qiime2.org/qiime2/2024.5/amplicon/released/qiime2-amplicon-2024.5-py38-linux-conda.yml"
YAML_FILE="qiime2-2024.5.yml"

# Detectar mamba o conda
if command -v mamba &> /dev/null; then
  CMD="mamba"
  echo ">> Usando mamba"
else
  CMD="conda"
  echo ">> Usando conda"
fi

echo ">> Eliminando entorno previo, si existe..."
$CMD env remove -n $ENV_NAME --yes >/dev/null 2>&1 || true

echo ">> Descargando YAML oficial de QIIME2 2024.5..."
curl -L $YAML_URL -o $YAML_FILE

if [ ! -s "$YAML_FILE" ]; then
  echo "ERROR: No se pudo descargar el YAML de QIIME2."
  exit 1
fi

echo ">> Creando entorno $ENV_NAME..."
$CMD env create -n $ENV_NAME --file $YAML_FILE

echo ">> Activando entorno..."
source $(dirname $(which $CMD))/../etc/profile.d/conda.sh
conda activate $ENV_NAME

echo ">> Instalando Deblur dentro del entorno..."
pip install --no-cache-dir deblur

echo ">> Verificando instalación..."
if ! command -v qiime &> /dev/null; then
  echo "ERROR: QIIME2 no se instaló correctamente."
  exit 1
fi

echo ">> Probando importación de Deblur en QIIME2..."
python - <<'EOF'
try:
    from qiime2.plugins.deblur.methods import denoise_16S
    print("OK: Deblur está disponible dentro de QIIME2.")
except Exception as e:
    print("ERROR: Deblur NO está disponible en QIIME2.")
    print(e)
    exit(1)
EOF

echo ""
echo "==========================================="
echo " QIIME2 2024.5 + Deblur instalado con éxito"
echo " Para usarlo:  conda activate $ENV_NAME"
echo "==========================================="
