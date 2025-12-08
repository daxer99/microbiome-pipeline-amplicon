#!/usr/bin/env bash
set -e

ENV_NAME="qiime2-amplicon-2024.2"
YAML_URL="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2024.2/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"

echo "======================================================"
echo " Instalación completa: QIIME2 2024.2 + Deblur + Plugins"
echo "======================================================"
echo ""

# -------------------------------------------------------
# 1. ELIMINAR MICROMAMBA DEL PATH
# -------------------------------------------------------
echo ">> Eliminando micromamba del PATH..."
export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v "micromamba" | paste -sd ":" -)

# -------------------------------------------------------
# 2. DETECTAR MAMBA O CONDA
# -------------------------------------------------------
if command -v mamba &>/dev/null; then
    CMD=mamba
    echo ">> Usando mamba"
else
    CMD=conda
    echo ">> Usando conda"
fi

# -------------------------------------------------------
# 3. FIX NECESARIO PARA QIIME2
# -------------------------------------------------------
echo ">> Ajustando channel_priority = flexible"
conda config --set channel_priority flexible

# -------------------------------------------------------
# 4. ELIMINAR ENTORNO PREVIO
# -------------------------------------------------------
echo ""
echo ">> Eliminando entorno anterior si existe..."
conda remove -n "$ENV_NAME" --all -y || true

# -------------------------------------------------------
# 5. DESCARGAR YAML (verificar)
# -------------------------------------------------------
echo ">> Descargando YAML oficial..."
curl -sSL "$YAML_URL" -o qiime2.yml

if ! grep -q "dependencies" qiime2.yml; then
    echo "ERROR: YAML inválido o no descargado correctamente."
    cat qiime2.yml
    exit 1
fi

# -------------------------------------------------------
# 6. CREAR ENTORNO SIN --solver
# -------------------------------------------------------
echo ">> Creando entorno QIIME2 (sin --solver)..."

$CMD env create -n "$ENV_NAME" -f qiime2.yml || {
    echo "Falló con mamba — reintentando con conda..."
    conda env create -n "$ENV_NAME" -f qiime2.yml
}

# -------------------------------------------------------
# 7. ACTIVAR ENTORNO
# -------------------------------------------------------
echo ">> Activando entorno..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# -------------------------------------------------------
# 8. INSTALAR DEBLUR + DEPENDENCIAS
# -------------------------------------------------------
echo ">> Instalando Deblur..."
pip install --upgrade pip
pip install sortmerna==2.0.0a
pip install deblur==1.1.1

# symlink requerido por Deblur
ln -sf "$(which sortmerna)" "$CONDA_PREFIX/bin/sortmerna"

# -------------------------------------------------------
# 9. PLUGINS EXTRA
# -------------------------------------------------------
echo ">> Instalando plugins extra..."
pip install \
    q2-sample-classifier \
    q2-quality-control \
    q2-fragment-insertion \
    q2-diversity-lib \
    q2-longitudinal

# -------------------------------------------------------
# 10. TESTS
# -------------------------------------------------------
echo ""
echo ">> Verificando instalación..."

qiime --help >/dev/null && echo "✓ qiime OK"
python -c "import deblur" && echo "✓ deblur OK"
qiime deblur denoise-16S --help >/dev/null && echo "✓ qiime deblur OK"

echo ""
echo "======================================================"
echo "  QIIME2 2024.2 + Deblur + Plugins INSTALADO 🎉"
echo "  Para activarlo: conda activate $ENV_NAME"
echo "======================================================"
