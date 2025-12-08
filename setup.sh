#!/usr/bin/env bash
set -e

echo "======================================"
echo "  Instalación QIIME2 Amplicon (latest) "
echo "======================================"

# Nombre del entorno
ENV_NAME="qiime2-amplicon-latest"

# URL del archivo .yml oficial para Linux
YML_URL="https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml"

YML_FILE="qiime2-amplicon-latest.yml"

echo "[1/5] Descargando definición del entorno..."
if command -v wget >/dev/null 2>&1; then
  wget -O ${YML_FILE} ${YML_URL} || { echo "ERROR: wget falló"; exit 1; }
elif command -v curl >/dev/null 2>&1; then
  curl -o ${YML_FILE} ${YML_URL} || { echo "ERROR: curl falló"; exit 1; }
else
  echo "ERROR: no tenés ni wget ni curl instalados. Instalalos primero."; exit 1
fi

echo "[OK] Archivo descargado: ${YML_FILE}"

echo "[2/5] Creando entorno conda \"${ENV_NAME}\"..."
conda env create -n ${ENV_NAME} --file ${YML_FILE}

echo "[3/5] Activando entorno..."
# Cargar configuración de shell si hace falta
if [[ $SHELL == *"bash"* ]]; then
  source ~/.bashrc 2>/dev/null || true
elif [[ $SHELL == *"zsh"* ]]; then
  source ~/.zshrc 2>/dev/null || true
fi
conda activate ${ENV_NAME}

echo "[4/5] Verificando instalación de QIIME2..."
qiime --help || { echo "ERROR: QIIME2 no se instaló correctamente"; exit 1; }

echo "[5/5] Instalación completada. Entorno activo: ${ENV_NAME}"
