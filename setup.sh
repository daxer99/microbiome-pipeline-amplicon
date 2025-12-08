#!/usr/bin/env bash
set -e

ENV_NAME="qiime2-amplicon-2024.2"
YAML_URL="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2024.2/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"

echo "======================================================"
echo " Instalación completa: QIIME2 2024.2 + Deblur + Plugins"
echo "======================================================"
echo ""

# -------------------------------------------------------
# 1. DESHABILITAR MICROMAMBA
# -------------------------------------------------------
echo ">> Eliminando micromamba del PATH temporalmente..."
export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v "micromamba" | paste -sd ":" -)

# -------------------------------------------------------
# 2. DETECTAR COMANDO (mamba > conda)
# -------------------------------------------------------
if command -v mamba &> /dev/null; then
    CMD=mamba
    echo ">> Usando mamba"
else
    CMD=conda
    echo ">> Usando conda"
fi

# -------------------------------------------------------
# 3. FIX DE QIIME2 (channel_priority flexible)
# -------------------------------------------------------
echo ""
echo ">> Ajustando channel_priority = flexible"
conda config --set channel_priority flexible

# -------------------------------------------------------
# 4. ELIMINAR ENTORNO PREVIO
# -------------------------------------------------------
echo ""
echo ">> Eliminando entorno previo si existe..."
conda remove -n "$ENV_NAME" --all -y || true

# -------------------------------------------------------
# 5. CREAR ENTORNO DESDE EL YAML OFICIAL DE QIIME2
# -------------------------------------------------------
echo ""
echo ">> Creando entorno con QIIME2 2024.2..."
$CMD env create \
    --name "$ENV_NAME" \
    --file "$YAML_URL" \
    --solver=libmamba

# -------------------------------------------------------
# 6. ACTIVAR ENTORNO
# -------------------------------------------------------
echo ""
echo ">> Activando entorno..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# -------------------------------------------------------
# 7. INSTALAR DEBLUR (pip) + DEPENDENCIAS
# -------------------------------------------------------
echo ""
echo ">> Instalando Deblur + dependencias..."
pip install --upgrade pip

# sortmerna 2 no existe en bioconda → parche compatible
pip install sortmerna==2.0.0a

# Deblur estable
pip install deblur==1.1.1

# FIX: crear symlink requerido por Deblur
mkdir -p $CONDA_PREFIX/bin
if [ ! -f $CONDA_PREFIX/bin/sortmerna ]; then
    ln -s $(which sortmerna) $CONDA_PREFIX/bin/sortmerna
fi

# -------------------------------------------------------
# 8. INSTALAR PLUGINS EXTRAS COMPATIBLES
# -------------------------------------------------------
echo ""
echo ">> Instalando plugins extra..."
pip install \
    q2-sample-classifier \
    q2-quality-control \
    q2-fragment-insertion \
    q2-diversity-lib \
    q2-longitudinal

# -------------------------------------------------------
# 9. TESTS DE INSTALACIÓN
# -------------------------------------------------------
echo ""
echo ">> Verificando QIIME2..."
qiime --help >/dev/null && echo "   ✓ qiime OK"

echo ">> Verificando Deblur..."
python -c "import deblur" && echo "   ✓ deblur importable"

echo ">> Verificando método: denoise-16S..."
qiime deblur denoise-16S --help >/dev/null && echo "   ✓ qiime deblur OK"

echo ""
echo "======================================================"
echo "  QIIME2 2024.2 INSTALADO CON DEBLUR + PLUGINS 🎉"
echo "  Para activarlo: conda activate $ENV_NAME"
echo "======================================================"
