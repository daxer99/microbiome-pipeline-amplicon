#!/bin/bash
# Microbiome Pipeline - Instalación Simplificada y Robusta
# Versión: 2.0

set -e  # Salir si hay error

echo "=============================================="
echo "   Microbiome Pipeline - Instalación v2.0"
echo "=============================================="
echo ""

# Colores para output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

info() { echo -e "${GREEN}[✓]${NC} $1"; }
warn() { echo -e "${YELLOW}[!]${NC} $1"; }
error() { echo -e "${RED}[✗]${NC} $1"; exit 1; }
step() { echo -e "${BLUE}[→]${NC} $1"; }

# ============================================
# 1. VERIFICACIONES PREVIAS
# ============================================

step "Verificando requisitos previos..."

# Verificar que estamos en el directorio correcto
if [ ! -f "environment.yml" ]; then
    error "No se encuentra environment.yml. Ejecuta este script desde el directorio del pipeline."
fi

# Verificar conda
if ! command -v conda &> /dev/null; then
    error "Conda no está instalado. Instala Miniconda desde: https://docs.conda.io/en/latest/miniconda.html"
fi
info "Conda encontrado: $(conda --version)"

# Verificar conexión a internet
if ! ping -c 1 8.8.8.8 &> /dev/null; then
    warn "No se detecta conexión a internet. La instalación podría fallar."
    read -p "¿Continuar de todos modos? (s/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Ss]$ ]]; then
        exit 1
    fi
fi

echo ""

# ============================================
# 2. LIMPIAR INSTALACIÓN PREVIA
# ============================================

step "Limpiando instalación previa (si existe)..."

# Desactivar cualquier entorno
conda deactivate 2>/dev/null || true

# Eliminar entorno anterior si existe
if conda env list | grep -q "qiime2-amplicon-2024.2"; then
    warn "Encontrado entorno previo. Eliminando..."
    conda env remove -n qiime2-amplicon-2024.2 -y
    info "Entorno anterior eliminado"
fi

# Limpiar caché de conda
conda clean --all -y &> /dev/null || true

echo ""

# ============================================
# 3. INSTALAR MAMBA (para instalación más rápida)
# ============================================

step "Verificando mamba..."

if ! command -v mamba &> /dev/null; then
    warn "Mamba no encontrado. Instalando para acelerar el proceso..."
    conda install -n base -c conda-forge mamba -y
    info "Mamba instalado"
else
    info "Mamba ya está instalado"
fi

echo ""

# ============================================
# 4. CREAR ENTORNO
# ============================================

step "Creando entorno (esto puede tomar 15-30 minutos)..."
echo ""
echo "   Instalando:"
echo "   • QIIME2 2024.2 con TODOS los plugins"
echo "   • PICRUSt2 para predicción funcional"
echo "   • Python 3.8 + librerías científicas"
echo "   • Herramientas bioinformáticas"
echo ""

# Usar mamba para instalación más rápida
if command -v mamba &> /dev/null; then
    info "Usando mamba para instalación acelerada..."
    if mamba env create -f environment.yml; then
        info "Entorno creado exitosamente con mamba"
    else
        warn "Mamba falló, intentando con conda..."
        conda env create -f environment.yml || error "Falló la creación del entorno"
    fi
else
    conda env create -f environment.yml || error "Falló la creación del entorno"
fi

echo ""

# ============================================
# 5. ACTIVAR ENTORNO
# ============================================

step "Activando entorno..."

# Inicializar conda para bash
eval "$(conda shell.bash hook)"

# Activar entorno
if conda activate qiime2-amplicon-2024.2; then
    info "Entorno activado correctamente"
else
    error "No se pudo activar el entorno"
fi

echo ""

# ============================================
# 6. INSTALAR SRA TOOLKIT (PARA DESCARGAS DE SRA)
# ============================================

step "Instalando SRA Toolkit para descarga de datos..."

# Intentar instalación con conda en el entorno activo
info "Instalando sra-tools desde bioconda..."
if conda install -c bioconda sra-tools -y; then
    info "SRA Toolkit instalado correctamente"
else
    warn "No se pudo instalar con bioconda, intentando alternativa..."

    # Intentar con mamba si está disponible
    if command -v mamba &> /dev/null; then
        mamba install -c bioconda sra-tools -y || warn "Falló instalación con mamba"
    fi
fi

# Verificar instalación
echo ""
step "Verificando instalación de SRA Toolkit..."
if command -v prefetch &> /dev/null; then
    info "prefetch disponible: $(prefetch --version 2>&1 | head -1 || echo 'instalado')"
else
    warn "prefetch no disponible - intentando instalación manual..."

    # Intentar instalar con pip como alternativa (si hay paquete pip)
    pip install sra-tools 2>/dev/null || true

    # Verificar de nuevo
    if command -v prefetch &> /dev/null; then
        info "prefetch ahora disponible"
    else
        warn "prefetch no instalado. Las funciones de descarga no funcionarán."
        warn "Para instalar manualmente después: conda install -c bioconda sra-tools"
    fi
fi

echo ""

# ============================================
# 7. INSTALAR DOKDO DESDE GITHUB (CORRECTO)
# ============================================

step "Instalando dokdo desde GitHub..."

info "Instalando dokdo desde el repositorio oficial..."
pip install git+https://github.com/sbslee/dokdo.git || warn "Falló instalación directa de dokdo"

# Verificar instalación
if python -c "import dokdo; print(f'✅ Dokdo importado exitosamente')" 2>/dev/null; then
    info "Dokdo instalado correctamente"
    # Obtener versión si es posible
    python -c "import dokdo; print(f'   Versión: {dokdo.__version__}')" 2>/dev/null || echo "   (versión no disponible)"
else
    warn "Dokdo no se pudo importar después de pip install"

    # Intentar método alternativo
    info "Intentando instalación desde clonado local..."
    DOKDO_TEMP_DIR=$(mktemp -d)
    cd "$DOKDO_TEMP_DIR"

    git clone https://github.com/sbslee/dokdo.git && cd dokdo && pip install .

    cd - > /dev/null
    rm -rf "$DOKDO_TEMP_DIR"

    # Verificar de nuevo
    if python -c "import dokdo" 2>/dev/null; then
        info "Dokdo instalado correctamente desde clonado"
    else
        warn "Dokdo no se pudo instalar. Algunas funciones de visualización no estarán disponibles."
        warn "Puedes instalarlo manualmente después con:"
        warn "  pip install git+https://github.com/sbslee/dokdo.git"
    fi
fi

echo ""

# ============================================
# 8. VERIFICAR INSTALACIÓN
# ============================================

step "Verificando instalación completa..."
echo ""

# Crear script de verificación temporal (INCLUYENDO DOKDO)
cat > /tmp/verify_installation.py << 'VERIFY_SCRIPT'
import sys
import subprocess
import os

print("Verificando componentes críticos:")
print("=" * 60)

# Lista de imports críticos que debe verificar
critical_imports = [
    ('qiime2', 'import qiime2'),
    ('Artifact', 'from qiime2 import Artifact, Metadata'),
    ('demux', 'from qiime2.plugins.demux.visualizers import summarize'),
    ('deblur', 'from qiime2.plugins.deblur.methods import denoise_16S'),
    ('feature-table', 'from qiime2.plugins.feature_table.methods import filter_features'),
    ('feature-classifier', 'from qiime2.plugins.feature_classifier.pipelines import classify_consensus_vsearch'),
    ('taxa', 'from qiime2.plugins.taxa.visualizers import barplot'),
    ('diversity', 'from qiime2.plugins.diversity.pipelines import alpha, beta'),
    ('phylogeny', 'from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree'),
    ('quality-control', 'from qiime2.plugins import quality_control'),
    ('pandas', 'import pandas'),
    ('click', 'import click'),
    ('biom', 'import biom'),
    # Dokdo - verificar pero no fallar si no está
    ('dokdo', 'import dokdo'),
]

failed_imports = []
for name, import_stmt in critical_imports:
    try:
        exec(import_stmt)
        print(f"✓ {name:20} OK")
    except Exception as e:
        if name == 'dokdo':
            print(f"⚠ {name:20} AUSENTE: {str(e)[:30]}")
        else:
            print(f"✗ {name:20} FALLO: {str(e)[:30]}")
            failed_imports.append(name)

print("")
print("Verificando herramientas de línea de comandos:")
print("-" * 40)

# Verificar herramientas de SRA
sra_tools = ['prefetch', 'fasterq-dump', 'vdb-validate']
missing_tools = []
for tool in sra_tools:
    try:
        # Verificar si el comando existe
        result = subprocess.run(['which', tool], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✓ {tool:15} OK")
        else:
            print(f"✗ {tool:15} NO ENCONTRADO")
            missing_tools.append(tool)
    except Exception:
        print(f"✗ {tool:15} ERROR AL VERIFICAR")
        missing_tools.append(tool)

print("=" * 60)

# Si dokdo está instalado, verificar funciones específicas
try:
    import dokdo
    dokdo_functions = ['quality_plot', 'taxa_barplot', 'alpha_diversity_plot']
    available_functions = []
    for func in dokdo_functions:
        if hasattr(dokdo, func):
            available_functions.append(func)

    if available_functions:
        print(f"✅ Dokdo funciones disponibles: {', '.join(available_functions)}")
    else:
        print("⚠  Dokdo instalado pero funciones no encontradas")
except ImportError:
    print("⚠  Dokdo no instalado - algunas visualizaciones no estarán disponibles")

print("=" * 60)

if not failed_imports:
    print("\n✅ COMPONENTES CRÍTICOS VERIFICADOS\n")
    sys.exit(0)
else:
    if failed_imports:
        print(f"\n⚠️  {len(failed_imports)} import(s) críticos fallaron: {', '.join(failed_imports)}")
    if missing_tools:
        print(f"⚠️  {len(missing_tools)} herramienta(s) faltan: {', '.join(missing_tools)}")
        print("   Para instalar: conda install -c bioconda sra-tools")
    print("")
    sys.exit(1)
VERIFY_SCRIPT

# Ejecutar verificación
if python /tmp/verify_installation.py; then
    info "Todos los componentes verificados correctamente"
else
    warn "Algunos componentes fallaron. Intentando reparar..."

    # Reinstalar paquetes problemáticos
    conda install -c https://packages.qiime2.org/qiime2/2024.2/amplicon/released \
        q2-feature-classifier q2-phylogeny q2-quality-control -y

    echo ""
    step "Verificando nuevamente después de reparación..."
    if python /tmp/verify_installation.py; then
        info "Reparación exitosa"
    else
        error "La instalación tiene problemas. Revisa los errores arriba."
    fi
fi

rm /tmp/verify_installation.py

echo ""

# ============================================
# 9. VERIFICAR PIPELINE CLI
# ============================================

step "Verificando pipeline CLI..."

if [ -f "microbiome_cli.py" ]; then
    if python microbiome_cli.py --help &> /dev/null; then
        info "Pipeline CLI funciona correctamente"
        echo ""
        echo "Comandos disponibles:"
        python microbiome_cli.py --help | grep -E "^  [a-z-]+" | head -10
    else
        warn "El CLI tiene problemas. Errores:"
        python microbiome_cli.py --help 2>&1 | head -5
    fi
else
    warn "microbiome_cli.py no encontrado en el directorio actual"
fi

echo ""

# ============================================
# 10. CREAR SCRIPTS DE AYUDA
# ============================================

step "Creando scripts de ayuda..."

# Script de activación
cat > activate.sh << 'EOF'
#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2

echo "=============================================="
echo "  ✅ Microbiome Pipeline Environment"
echo "=============================================="
echo ""
echo "Para usar el pipeline:"
echo "  python microbiome_cli.py --help"
echo ""
echo "Ejemplo rápido:"
echo "  python microbiome_cli.py download samples.csv"
echo ""
EOF
chmod +x activate.sh

# Script de test rápido (actualizado)
cat > test.sh << 'EOF'
#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.2

echo "=============================================="
echo "  🧪 Test Rápido - Microbiome Pipeline"
echo "=============================================="
echo ""
echo "1. Probando imports críticos..."
python -c "
from qiime2 import Artifact
from qiime2.plugins.demux.visualizers import summarize
from qiime2.plugins.feature_classifier.pipelines import classify_consensus_vsearch
import pandas
try:
    import dokdo
    print('   ✅ Dokdo instalado')
except:
    print('   ⚠  Dokdo no instalado (algunas visualizaciones no funcionarán)')
print('   ✅ Todos los imports funcionan')
"

echo ""
echo "2. Probando herramientas SRA..."
for tool in prefetch fasterq-dump vdb-validate; do
    if command -v $tool &> /dev/null; then
        echo "   ✅ $tool disponible"
    else
        echo "   ❌ $tool NO disponible"
    fi
done

echo ""
echo "3. Probando CLI..."
python microbiome_cli.py --help | head -3
echo ""
echo "=============================================="
echo ""
EOF
chmod +x test.sh

info "Scripts creados: activate.sh y test.sh"

echo ""

# ============================================
# 11. RESUMEN FINAL
# ============================================

echo "=============================================="
echo "  ✅ INSTALACIÓN COMPLETADA"
echo "=============================================="
echo ""
echo "Para usar el pipeline:"
echo ""
echo "  1. Activar entorno:"
echo "     source activate.sh"
echo "     # o manualmente:"
echo "     conda activate qiime2-amplicon-2024.2"
echo ""
echo "  2. Ejecutar pipeline:"
echo "     python microbiome_cli.py --help"
echo ""
echo "  3. Test rápido:"
echo "     ./test.sh"
echo ""
echo "  4. Descargar datos de ejemplo:"
echo "     python microbiome_cli.py download example/sra.csv --output-dir example/"
echo ""
echo "=============================================="
echo ""

# Información de versiones instaladas
step "Versiones instaladas:"
echo "  • QIIME2: $(python -c 'import qiime2; print(qiime2.__version__)' 2>/dev/null || echo 'error')"
echo "  • Python: $(python --version 2>&1 | cut -d' ' -f2)"
echo "  • PICRUSt2: $(picrust2_pipeline.py -v 2>&1 | head -1 || echo 'instalado')"
echo "  • SRA Toolkit: $(prefetch --version 2>&1 | grep -o 'version [0-9.]*' | head -1 || echo 'instalado')"
echo "  • Dokdo: $(python -c 'import dokdo; print(dokdo.__version__)' 2>/dev/null || echo 'no instalado')"

echo ""
info "Todo listo para usar el pipeline"
echo ""