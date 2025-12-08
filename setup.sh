#!/bin/bash
# Script de instalación limpia para Microbiome Pipeline

echo "=================================================="
echo "   INSTALACIÓN LIMPIA - Microbiome Pipeline"
echo "=================================================="
echo ""

# Configuración
ENV_NAME="qiime2-amplicon-2024.2"
PYTHON_VERSION="3.8"

# Colores para mensajes
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# ===================================================
# 1. VERIFICAR MINIFORGE3
# ===================================================
echo -e "${BLUE}[1] Verificando Miniforge3...${NC}"

# Verificar que estamos usando miniforge3
if [[ "$(which conda)" != *"miniforge3"* ]]; then
    echo -e "${RED}ERROR: No se detecta Miniforge3 como el conda activo.${NC}"
    echo "Para activar Miniforge3, ejecuta:"
    echo "  source ~/miniforge3/etc/profile.d/conda.sh"
    echo "O añade esta línea a tu ~/.bashrc:"
    echo "  source ~/miniforge3/etc/profile.d/conda.sh"
    exit 1
fi

echo "✓ Usando Miniforge3 en: $(which conda)"
echo ""

# ===================================================
# 2. CONFIGURAR CANALES
# ===================================================
echo -e "${BLUE}[2] Configurando canales de conda...${NC}"

# Configurar canales en el orden correcto
conda config --remove-key channels 2>/dev/null || true
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels qiime2
conda config --set channel_priority strict

echo "✓ Canales configurados: defaults, bioconda, conda-forge, qiime2"
echo ""

# ===================================================
# 3. INSTALAR/ACTUALIZAR MAMBA
# ===================================================
echo -e "${BLUE}[3] Instalando Mamba...${NC}"

# Verificar si mamba está instalado
if ! command -v mamba &> /dev/null; then
    echo "Instalando mamba..."
    conda install -n base -c conda-forge mamba -y
fi

# Actualizar mamba si ya está instalado
mamba update -n base -c conda-forge mamba -y

echo "✓ Mamba instalado/actualizado: $(mamba --version)"
echo ""

# ===================================================
# 4. CREAR ENTORNO NUEVO
# ===================================================
echo -e "${BLUE}[4] Creando entorno '$ENV_NAME'...${NC}"

# Primero, verificar si el entorno ya existe
if conda env list | grep -q "^$ENV_NAME\s"; then
    echo "⚠  El entorno ya existe. Eliminando..."
    conda env remove -n "$ENV_NAME" -y
fi

# Crear entorno con mamba (más rápido)
echo "Creando entorno con Python $PYTHON_VERSION y QIIME2 2024.2..."
echo "Esto puede tomar 15-30 minutos..."

mamba create -n "$ENV_NAME" \
    python="$PYTHON_VERSION" \
    qiime2-amplicon=2024.2 \
    snakemake-minimal=7.32.4 \
    multiqc=1.17 \
    fastp=0.23.4 \
    pandas=2.0.3 \
    numpy=1.24.3 \
    matplotlib=3.7.2 \
    seaborn=0.12.2 \
    sra-tools \
    -c qiime2 \
    -c conda-forge \
    -c bioconda \
    -c defaults \
    -y

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Entorno creado exitosamente${NC}"
else
    echo -e "${RED}✗ Error al crear el entorno${NC}"
    exit 1
fi

echo ""

# ===================================================
# 5. ACTIVAR ENTORNO E INSTALAR DOKDO
# ===================================================
echo -e "${BLUE}[5] Activando entorno e instalando dependencias...${NC}"

# Activar el entorno
eval "$(conda shell.bash hook)"
mamba activate "$ENV_NAME"

if [[ "$CONDA_DEFAULT_ENV" != "$ENV_NAME" ]]; then
    echo "⚠  No se pudo activar automáticamente. Intentando manualmente..."
    conda activate "$ENV_NAME"
fi

echo "✓ Entorno activado: $(python --version)"

# Instalar dokdo
echo "Instalando dokdo..."
pip install --no-deps git+https://github.com/sbslee/dokdo.git 2>/dev/null || \
    pip install git+https://github.com/sbslee/dokdo.git 2>/dev/null || \
    echo "⚠  No se pudo instalar dokdo (continuando sin él)"

echo ""

# ===================================================
# 6. VERIFICAR INSTALACIÓN
# ===================================================
echo -e "${BLUE}[6] Verificando instalación...${NC}"

# Función de verificación
verify_tool() {
    local tool=$1
    local command=$2
    if command -v $tool &> /dev/null; then
        echo -e "  ${GREEN}✓${NC} $tool"
        return 0
    else
        echo -e "  ${RED}✗${NC} $tool (no encontrado)"
        return 1
    fi
}

echo "Verificando herramientas:"
verify_tool "qiime" "qiime --version"
verify_tool "snakemake" "snakemake --version"
verify_tool "fastp" "fastp --version"
verify_tool "multiqc" "multiqc --version"
verify_tool "python" "python --version"

# Verificar QIIME2 específicamente
echo ""
echo "Verificando QIIME2:"
if python -c "import qiime2; print(f'  ✓ QIIME2 {qiime2.__version__}')" 2>/dev/null; then
    echo -e "  ${GREEN}✓ QIIME2 cargado correctamente${NC}"
else
    echo -e "  ${RED}✗ No se pudo importar QIIME2${NC}"
fi

echo ""

# ===================================================
# 7. CONFIGURAR PIPELINE
# ===================================================
echo -e "${BLUE}[7] Configurando estructura del pipeline...${NC}"

# Crear directorios necesarios
mkdir -p data/raw data/processed logs reports config scripts

# Copiar configuración si existe
if [ -f "config/config_template.yaml" ] && [ ! -f "config/config.yaml" ]; then
    cp config/config_template.yaml config/config.yaml
    echo "✓ Archivo de configuración creado: config/config.yaml"
else
    echo "⚠  No se encontró config/config_template.yaml"
fi

# Hacer ejecutables los scripts
chmod +x *.py *.sh 2>/dev/null || true

echo "✓ Estructura de directorios creada"
echo ""

# ===================================================
# 8. CREAR SCRIPTS DE AYUDA
# ===================================================
echo -e "${BLUE}[8] Creando scripts de ayuda...${NC}"

# Script de activación
cat > activate_microbiome.sh << 'EOF'
#!/bin/bash
# Script para activar el entorno del pipeline

# Inicializar conda
source ~/miniforge3/etc/profile.d/conda.sh

# Activar entorno
mamba activate qiime2-amplicon-2024.2

# Verificar activación
if [[ "$CONDA_DEFAULT_ENV" == "qiime2-amplicon-2024.2" ]]; then
    echo "=============================================="
    echo "  ✅ Entorno qiime2-amplicon-2024.2 activado"
    echo "=============================================="
    echo ""
    echo "Herramientas disponibles:"
    echo "  • qiime --help"
    echo "  • snakemake --help"
    echo "  • fastp --help"
    echo "  • multiqc --help"
    echo ""
    echo "Para ejecutar el pipeline:"
    echo "  python microbiome_cli.py --help"
    echo ""
else
    echo "❌ No se pudo activar el entorno"
    exit 1
fi
EOF

chmod +x activate_microbiome.sh

# Script de prueba
cat > test_install.sh << 'EOF'
#!/bin/bash
# Script de prueba de instalación

echo "🧪 Probando instalación del Microbiome Pipeline"
echo "=============================================="

# Activar entorno
source activate_microbiome.sh

echo ""
echo "1. Verificando versiones:"
echo "   Python: $(python --version 2>&1)"
echo "   QIIME2: $(python -c 'import qiime2; print(qiime2.__version__)' 2>/dev/null || echo 'No disponible')"
echo "   Snakemake: $(snakemake --version 2>/dev/null | head -1 || echo 'No disponible')"

echo ""
echo "2. Verificando directorios:"
for dir in data/raw data/processed logs reports config; do
    if [ -d "$dir" ]; then
        echo "   ✅ $dir"
    else
        echo "   ❌ $dir (faltante)"
    fi
done

echo ""
echo "3. Probar importación de módulos:"
python -c "
import sys
modules = ['qiime2', 'pandas', 'numpy', 'seaborn']
for mod in modules:
    try:
        __import__(mod)
        print(f'   ✅ {mod}')
    except:
        print(f'   ❌ {mod}')
"

echo ""
echo "✅ Prueba completada"
EOF

chmod +x test_install.sh

echo "✓ Scripts creados: activate_microbiome.sh, test_install.sh"
echo ""

# ===================================================
# 9. RESUMEN FINAL
# ===================================================
echo "=================================================="
echo "   🎉 INSTALACIÓN COMPLETADA"
echo "=================================================="
echo ""
echo "RESUMEN:"
echo "  • Entorno: qiime2-amplicon-2024.2"
echo "  • Python: $PYTHON_VERSION"
echo "  • Gestor: Miniforge3 + Mamba"
echo "  • Ubicación: ~/miniforge3/envs/qiime2-amplicon-2024.2"
echo ""
echo "PARA USAR EL PIPELINE:"
echo ""
echo "  1. Activar el entorno:"
echo "     source activate_microbiome.sh"
echo ""
echo "  2. Probar la instalación:"
echo "     ./test_install.sh"
echo ""
echo "  3. Configurar el pipeline (edita):"
echo "     config/config.yaml"
echo ""
echo "  4. Ejecutar el pipeline:"
echo "     python microbiome_cli.py --help"
echo ""
echo "  5. Para análisis interactivo:"
echo "     jupyter lab"
echo ""
echo "=================================================="
echo ""
echo "📌 NOTAS:"
echo "  • El entorno está aislado y no interfiere con otras instalaciones"
echo "  • Para desactivar: conda deactivate"
echo "  • Para eliminar: conda env remove -n qiime2-amplicon-2024.2"
echo ""