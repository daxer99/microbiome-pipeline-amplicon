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

echo ">> Eliminando entorno previo si existe..."
$CMD env remove -n $ENV_NAME --yes >/dev/null 2>&1 || true

echo ">> Descargando YAML oficial..."
curl -L $YAML_URL -o $YAML_FILE

echo ">> Verificando que el YAML sea válido..."
if ! grep -q "dependencies:" "$YAML_FILE"; then
    echo "ERROR: Archivo YAML inválido o corrupto."
    echo "Descargado:"
    head "$YAML_FILE"
    exit 1
fi

echo ">> Creando entorno $ENV_NAME..."
$CMD env create -n $ENV_NAME --file $YAML_FILE

echo ">> Activando entorno..."
source $(dirname $(which conda))/../etc/profile.d/conda.sh
conda activate $ENV_NAME

echo ">> Instalando Deblur..."
pip install --no-cache-dir deblur

echo ">> Probando import de Deblur en QIIME2..."
python - <<'EOF'
try:
    from qiime2.plugins.deblur.methods import denoise_16S
    print("OK: Deblur está disponible.")
except Exception as e:
    print("ERROR: Deblur no se cargó.")
    print(e)
    exit(1)
EOF

echo "==========================================="
echo " QIIME2 2024.5 + Deblur instalado con éxito"
echo " Activar con: conda activate $ENV_NAME"
echo "==========================================="
