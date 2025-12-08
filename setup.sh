#!/usr/bin/env bash
set -e

echo "======================================"
echo "  Instalación QIIME2 2024.2 OFICIAL"
echo "======================================"

# 1. Descargar archivo oficial de entorno
echo "[1/4] Descargando entorno oficial..."
wget -q https://data.qiime2.org/distro/core/qiime2-amplicon-2024.2-py38-linux-conda.yml

echo "[OK] Archivo descargado: qiime2-amplicon-2024.2-py38-linux-conda.yml"

# 2. Crear entorno
echo "======================================"
echo "[2/4] Creando entorno QIIME2..."
echo "======================================"

conda env create -n qiime2-2024.2 --file qiime2-amplicon-2024.2-py38-linux-conda.yml

echo "[OK] Entorno qiime2-2024.2 creado."

# 3. Activar entorno
echo "======================================"
echo "[3/4] Activando entorno..."
echo "======================================"

# Cargar configuración del shell
if [[ $SHELL == *"bash"* ]]; then
    source ~/.bashrc 2>/dev/null || true
elif [[ $SHELL == *"zsh"* ]]; then
    source ~/.zshrc 2>/dev/null || true
fi

conda activate qiime2-2024.2

echo "[OK] Entorno activado."

# 4. Verificar funcionamiento
echo "======================================"
echo "[4/4] Verificando QIIME2"
echo "======================================"

qiime --help || { echo "ERROR: QIIME2 no se instaló correctamente"; exit 1; }

echo "======================================"
echo "  QIIME2 2024.2 instalado correctamente"
echo "======================================"
echo "  Entorno activo: qiime2-2024.2"
echo "======================================"