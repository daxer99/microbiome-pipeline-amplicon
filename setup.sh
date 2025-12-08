#!/usr/bin/env bash
set -e

echo "===================================="
echo " CONFIGURANDO QIIME2 2024.10 + DEBLUR"
echo "===================================="

# -------------------------------------------------------------------
# Detectar gestor: mamba > conda
# -------------------------------------------------------------------
if command -v mamba &> /dev/null; then
    CMD="mamba"
    echo ">> Gestor detectado: mamba"
elif command -v conda &> /dev/null; then
    CMD="conda"
    echo ">> Gestor detectado: conda"
else
    echo "ERROR: No se encontró conda ni mamba en el sistema."
    exit 1
fi

# -------------------------------------------------------------------
# Nombre del entorno
# -------------------------------------------------------------------
ENV_NAME="qiime2-amplicon-2024.10"

# -------------------------------------------------------------------
# URL oficial del entorno QIIME2 Amplicon 2024.10
# -------------------------------------------------------------------
YAML_URL="https://raw.githubusercontent.com/qiime2/distributions/2024.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"

echo ">> Creando entorno $ENV_NAME usando YAML oficial..."
$CMD env create \
    --name "$ENV_NAME" \
    --file "$YAML_URL" || {
        echo "El entorno ya existe. Continuando..."
    }

# -------------------------------------------------------------------
# Activar entorno correctamente
# -------------------------------------------------------------------
echo ">> Activando entorno..."
if [ "$CMD" = "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
else
    eval "$(conda shell.bash hook)"
fi

$CMD activate "$ENV_NAME"

echo ">> Entorno activado: $ENV_NAME"

# -------------------------------------------------------------------
# Instalación de Deblur (bioconda)
# Nota: Amplicon ya incluye Deblur, pero forzamos instalación adicional por compatibilidad.
# -------------------------------------------------------------------
echo ">> Instalando Deblur desde bioconda..."
$CMD install -y -c bioconda -c conda-forge deblur

echo ""
echo "===================================="
echo "  INSTALACIÓN COMPLETA"
echo "===================================="
echo "Para activar QIIME2 en nuevas sesiones:"
echo ""
echo "    $CMD activate $ENV_NAME"
echo ""
