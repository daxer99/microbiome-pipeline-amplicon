#!/bin/bash

# Colores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Nombre del entorno
ENV_NAME="qiime2-amplicon-2024.10"

echo "=============================================="
echo "   Instalación de Microbiome Pipeline"
echo "   QIIME2 Amplicon 2024.10"
echo "=============================================="
echo ""

# Verificar que conda está instalado
if ! command -v conda &> /dev/null; then
    echo -e "${RED}❌ Error: conda no está instalado${NC}"
    echo "Por favor instala Miniconda o Anaconda primero:"
    echo "https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo -e "${GREEN}✓${NC} conda encontrado"

# Actualizar conda
echo ""
echo "Actualizando conda..."
conda update -n base -c defaults conda -y

# Eliminar entorno anterior si existe
if conda env list | grep -q "^${ENV_NAME} "; then
    echo ""
    echo -e "${YELLOW}⚠ Entorno '${ENV_NAME}' ya existe${NC}"
    read -p "¿Deseas eliminarlo y crear uno nuevo? (s/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Ss]$ ]]; then
        echo "Eliminando entorno anterior..."
        conda env remove -n ${ENV_NAME} -y
    else
        echo "Cancelando instalación."
        exit 0
    fi
fi

# Detectar sistema operativo
OS_TYPE=$(uname -s)
ARCH_TYPE=$(uname -m)

echo ""
echo "Sistema detectado: ${OS_TYPE} ${ARCH_TYPE}"

# Determinar URL del archivo de entorno según el sistema
if [[ "$OS_TYPE" == "Linux" ]]; then
    QIIME_URL="https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    echo "Usando configuración para Linux"
elif [[ "$OS_TYPE" == "Darwin" ]]; then
    if [[ "$ARCH_TYPE" == "arm64" ]]; then
        # Mac con Apple Silicon (M1/M2/M3)
        echo -e "${YELLOW}Detectado Mac con Apple Silicon${NC}"
        echo "Configurando instalación en modo Rosetta 2..."
        CONDA_SUBDIR=osx-64
        export CONDA_SUBDIR
        QIIME_URL="https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml"
    else
        # Mac Intel
        QIIME_URL="https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml"
        echo "Usando configuración para Mac Intel"
    fi
else
    echo -e "${RED}❌ Sistema operativo no soportado: ${OS_TYPE}${NC}"
    exit 1
fi

# Crear entorno QIIME2
echo ""
echo "=============================================="
echo "Creando entorno ${ENV_NAME}..."
echo "Esto puede tardar 10-20 minutos..."
echo "=============================================="
echo ""

if conda env create -n ${ENV_NAME} --file ${QIIME_URL}; then
    echo -e "${GREEN}✓ Entorno QIIME2 creado exitosamente${NC}"
else
    echo -e "${RED}❌ Error al crear el entorno QIIME2${NC}"
    exit 1
fi

# Activar entorno
echo ""
echo "Activando entorno..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${ENV_NAME}

if [[ $CONDA_DEFAULT_ENV != ${ENV_NAME} ]]; then
    echo -e "${RED}❌ Error al activar el entorno${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Entorno activado${NC}"

# Instalar dependencias adicionales
echo ""
echo "=============================================="
echo "Instalando dependencias adicionales..."
echo "=============================================="
echo ""

# PICRUSt2
echo "Instalando PICRUSt2..."
if conda install -c bioconda -c conda-forge picrust2 -y; then
    echo -e "${GREEN}✓ PICRUSt2 instalado${NC}"
else
    echo -e "${YELLOW}⚠ Advertencia: Error al instalar PICRUSt2${NC}"
fi

# SRA Toolkit para descarga de datos
echo ""
echo "Instalando SRA Toolkit..."
if conda install -c bioconda sra-tools -y; then
    echo -e "${GREEN}✓ SRA Toolkit instalado${NC}"
else
    echo -e "${YELLOW}⚠ Advertencia: Error al instalar SRA Toolkit${NC}"
fi

# Paquetes de Python adicionales
echo ""
echo "Instalando paquetes de Python adicionales..."
if pip install dokdo==1.16.0; then
    echo -e "${GREEN}✓ dokdo instalado${NC}"
else
    echo -e "${YELLOW}⚠ Advertencia: Error al instalar dokdo${NC}"
fi

# Verificar instalación
echo ""
echo "=============================================="
echo "Verificando instalación..."
echo "=============================================="
echo ""

# Verificar QIIME2
if qiime --version &> /dev/null; then
    QIIME_VERSION=$(qiime --version 2>&1)
    echo -e "${GREEN}✓ QIIME2: ${QIIME_VERSION}${NC}"
else
    echo -e "${RED}❌ Error: QIIME2 no se instaló correctamente${NC}"
    exit 1
fi

# Verificar PICRUSt2
if command -v picrust2_pipeline.py &> /dev/null; then
    PICRUST_VERSION=$(picrust2_pipeline.py --version 2>&1 | head -n 1)
    echo -e "${GREEN}✓ PICRUSt2: ${PICRUST_VERSION}${NC}"
else
    echo -e "${YELLOW}⚠ PICRUSt2 no disponible${NC}"
fi

# Verificar fastq-dump
if command -v fastq-dump &> /dev/null; then
    echo -e "${GREEN}✓ SRA Toolkit instalado${NC}"
else
    echo -e "${YELLOW}⚠ SRA Toolkit no disponible${NC}"
fi

# Verificar el CLI de la aplicación
if [[ -f "microbiome_cli.py" ]]; then
    if python microbiome_cli.py --help &> /dev/null; then
        echo -e "${GREEN}✓ microbiome_cli.py funcional${NC}"
    else
        echo -e "${YELLOW}⚠ Hay problemas con microbiome_cli.py${NC}"
    fi
else
    echo -e "${YELLOW}⚠ microbiome_cli.py no encontrado${NC}"
fi

# Instrucciones finales
echo ""
echo "=============================================="
echo -e "${GREEN}✓ ¡Instalación completada!${NC}"
echo "=============================================="
echo ""
echo "Para usar el pipeline:"
echo ""
echo "1. Activa el entorno:"
echo -e "   ${GREEN}conda activate ${ENV_NAME}${NC}"
echo ""
echo "2. Verifica la instalación:"
echo -e "   ${GREEN}qiime --help${NC}"
echo -e "   ${GREEN}python microbiome_cli.py --help${NC}"
echo ""
echo "3. Para desactivar el entorno:"
echo -e "   ${GREEN}conda deactivate${NC}"
echo ""
echo "4. Para eliminar el entorno (si es necesario):"
echo -e "   ${GREEN}conda env remove -n ${ENV_NAME}${NC}"
echo ""
echo "=============================================="
echo "Documentación y recursos:"
echo "  - QIIME2: https://docs.qiime2.org/2024.10/"
echo "  - Tu repositorio: https://github.com/daxer99/microbiome-pipeline-amplicon"
echo "=============================================="