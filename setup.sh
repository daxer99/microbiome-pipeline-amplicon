#!/usr/bin/env bash
set -e

echo "======================================================"
echo " Instalación completa: QIIME2 2024.2 + Deblur + Plugins"
echo "======================================================"

# -------------------------------------------------------------------
# 1) Usar SIEMPRE el mamba correcto (el de miniforge, NO micromamba)
# -------------------------------------------------------------------
MAMBA="/home/rodrigo/miniforge3/bin/mamba"
CONDA="/home/rodrigo/miniforge3/bin/conda"

echo ">> Usando mamba real: $MAMBA"

# -------------------------------------------------------------------
# 2) channel_priority flexible (soluciona el problema oficial del foro)
# -------------------------------------------------------------------
echo ">> Ajustando channel_priority = flexible"
$CONDA config --set channel_priority flexible

# -------------------------------------------------------------------
# 3) Borrar entorno previo dentro de Miniforge
# -------------------------------------------------------------------
echo
echo ">> Eliminando entorno previo si existe (solo en Miniforge)..."

$CONDA remove -n qiime2-amplicon-2024.2 --all -y || true


# -------------------------------------------------------------------
# 4) Crear entorno QIIME2 2024.2 desde el YAML oficial
# -------------------------------------------------------------------
echo
echo ">> Creando entorno con QIIME2 2024.2 en Miniforge..."

$MAMBA env create \
    -n qiime2-amplicon-2024.2 \
    --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2024.2/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml


# -------------------------------------------------------------------
# 5) Activar entorno correctamente
# -------------------------------------------------------------------
echo
echo ">> Activando entorno para instalar plugins..."

source /home/rodrigo/miniforge3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.2


# -------------------------------------------------------------------
# 6) Instalar Deblur (vía pip, como recomienda QIIME2)
# -------------------------------------------------------------------
echo
echo ">> Instalando Deblur..."
pip install deblur


# -------------------------------------------------------------------
# 7) Instalar plugins adicionales
# -------------------------------------------------------------------
echo
echo ">> Instalando plugins extra..."


echo " - q2-picrust2"
pip install git+https://github.com/picrust/picrust2.git
pip install git+https://github.com/biocore/q2-picrust2.git

echo " - q2-sourmash"
pip install git+https://github.com/sourmash-bio/q2-sourmash.git

echo " - q2-kaiju"
pip install git+https://github.com/bokulich-lab/q2-kaiju.git


# -------------------------------------------------------------------
# 8) Verificación final
# -------------------------------------------------------------------
echo
echo ">> Verificando instalación final..."
which qiime
qiime --help >/dev/null && echo "✔ QIIME2 funciona correctamente"


echo
echo "======================================================"
echo " Instalación finalizada CON ÉXITO"
echo " Para activar el entorno:"
echo "     conda activate qiime2-amplicon-2024.2"
echo "======================================================"
