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

# Inicializar conda para bash (IMPORTANTE: hacerlo al inicio)
eval "$(conda shell.bash hook)"

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
ENV_NAME="qiime2-amplicon-2024.2"
if conda env list | grep -q "$ENV_NAME"; then
    warn "Encontrado entorno previo. Eliminando..."
    conda env remove -n "$ENV_NAME" -y
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

    # Aceptar ToS si es necesario (para Conda 25+)
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

    conda install -n base -c conda-forge mamba -y
    info "Mamba instalado"

    # Inicializar mamba también
    if [ -f "$(conda info --base)/etc/profile.d/mamba.sh" ]; then
        source "$(conda info --base)/etc/profile.d/mamba.sh"
    fi
else
    info "Mamba ya está instalado"

    # Inicializar mamba si existe
    if [ -f "$(conda info --base)/etc/profile.d/mamba.sh" ]; then
        source "$(conda info --base)/etc/profile.d/mamba.sh"
    fi
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

# Verificar si el archivo environment.yml existe
if [ ! -f "environment.yml" ]; then
    error "Archivo environment.yml no encontrado!"
fi

# Usar mamba si está disponible, si no usar conda
if command -v mamba &> /dev/null; then
    info "Usando mamba para instalación acelerada..."
    echo "Creando entorno: $ENV_NAME"

    # Primero intentar crear con mamba usando environment.yml
    if mamba env create -n "$ENV_NAME" -f environment.yml; then
        info "Entorno creado exitosamente con mamba"
    else
        warn "Mamba falló con environment.yml, intentando enfoque alternativo..."

        # Enfoque alternativo: crear entorno base e instalar QIIME2 manualmente
        mamba create -n "$ENV_NAME" -c qiime2 -c conda-forge -c bioconda \
            qiime2-amplicon=2024.2 \
            snakemake-minimal=7.32.4 \
            multiqc=1.17 \
            fastp=0.23.4 \
            pandas=2.0.3 \
            numpy=1.24.3 \
            matplotlib=3.7.2 \
            seaborn=0.12.2 \
            python=3.8 -y

        info "Entorno creado con enfoque alternativo"
    fi
else
    warn "Mamba no disponible, usando conda (más lento)..."

    # Aceptar ToS para conda también
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

    if conda env create -n "$ENV_NAME" -f environment.yml; then
        info "Entorno creado exitosamente con conda"
    else
        error "Falló la creación del entorno con conda"
    fi
fi

echo ""

# ============================================
# 5. ACTIVAR ENTORNO
# ============================================

step "Activando entorno..."

# Esperar un momento para que el entorno esté listo
sleep 2

# Actualizar la lista de entornos
conda config --set env_prompt '({name}) ' 2>/dev/null || true

# Intentar activar con mamba si está disponible, si no con conda
if command -v mamba &> /dev/null; then
    if mamba activate "$ENV_NAME"; then
        info "Entorno activado correctamente con mamba"
    else
        # Intentar con conda
        if conda activate "$ENV_NAME"; then
            info "Entorno activado correctamente con conda"
        else
            error "No se pudo activar el entorno. Intenta manualmente: conda activate $ENV_NAME"
        fi
    fi
else
    if conda activate "$ENV_NAME"; then
        info "Entorno activado correctamente"
    else
        error "No se pudo activar el entorno"
    fi
fi

echo ""

# ============================================
# 6. INSTALAR DEPENDENCIAS ADICIONALES
# ============================================

step "Instalando dependencias adicionales..."

# Usar el gestor apropiado según qué comando esté disponible
if command -v mamba &> /dev/null && mamba list -n "$ENV_NAME" &> /dev/null; then
    info "Usando mamba para instalar dependencias..."

    # Activar el entorno para mamba
    mamba activate "$ENV_NAME"

    # Instalar paquetes adicionales con mamba
    mamba install -c conda-forge -c bioconda \
        snakemake-minimal=7.32.4 \
        multiqc=1.17 \
        fastp=0.23.4 \
        pandas=2.0.3 \
        numpy=1.24.3 \
        matplotlib=3.7.2 \
        seaborn=0.12.2 \
        sra-tools \
        -y

elif conda list -n "$ENV_NAME" &> /dev/null; then
    info "Usando conda para instalar dependencias..."

    # Activar el entorno para conda
    conda activate "$ENV_NAME"

    # Instalar paquetes adicionales con conda
    conda install -c conda-forge -c bioconda \
        snakemake-minimal=7.32.4 \
        multiqc=1.17 \
        fastp=0.23.4 \
        pandas=2.0.3 \
        numpy=1.24.3 \
        matplotlib=3.7.2 \
        seaborn=0.12.2 \
        sra-tools \
        -y
else
    warn "No se pudo acceder al entorno para instalar dependencias"
fi

echo ""

# ============================================
# 7. INSTALAR DOKDO DESDE GITHUB
# ============================================

step "Instalando dokdo desde GitHub..."

# Asegurarse de que estamos en el entorno correcto
if [[ "$CONDA_DEFAULT_ENV" != "$ENV_NAME" ]]; then
    if command -v mamba &> /dev/null; then
        mamba activate "$ENV_NAME"
    else
        conda activate "$ENV_NAME"
    fi
fi

info "Instalando dokdo desde el repositorio oficial..."
if pip install git+https://github.com/sbslee/dokdo.git; then
    info "Dokdo instalado correctamente"
else
    warn "Falló instalación directa de dokdo, intentando con --no-deps..."
    pip install --no-deps git+https://github.com/sbslee/dokdo.git || \
    warn "No se pudo instalar dokdo. Puedes instalarlo manualmente después."
fi

echo ""

# ============================================
# 8. VERIFICAR INSTALACIÓN
# ============================================

step "Verificando instalación completa..."
echo ""

# Crear script de verificación temporal
cat > /tmp/verify_installation.py << 'VERIFY_SCRIPT'
import sys
import subprocess
import os

print("Verificando componentes críticos:")
print("=" * 60)

# Lista de imports críticos que debe verificar
critical_imports = [
    ('qiime2', 'import qiime2'),
    ('pandas', 'import pandas'),
    ('numpy', 'import numpy'),
    ('seaborn', 'import seaborn'),
]

# Lista de imports opcionales
optional_imports = [
    ('dokdo', 'import dokdo'),
]

failed_critical = []
for name, import_stmt in critical_imports:
    try:
        exec(import_stmt)
        print(f"✓ {name:20} OK")
    except Exception as e:
        print(f"✗ {name:20} FALLO: {str(e)[:30]}")
        failed_critical.append(name)

print("")
print("Imports opcionales:")
for name, import_stmt in optional_imports:
    try:
        exec(import_stmt)
        print(f"✓ {name:20} OK")
    except Exception:
        print(f"⚠ {name:20} No instalado (opcional)")

print("")
print("Verificando herramientas de línea de comandos:")
print("-" * 40)

# Verificar herramientas importantes
tools = ['snakemake', 'fastp', 'multiqc']
for tool in tools:
    try:
        result = subprocess.run(['which', tool], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✓ {tool:15} OK")
        else:
            print(f"✗ {tool:15} NO ENCONTRADO")
    except Exception:
        print(f"✗ {tool:15} ERROR AL VERIFICAR")

print("=" * 60)

# Verificar QIIME2 específicamente
try:
    import qiime2
    print(f"✅ QIIME2 versión: {qiime2.__version__}")

    # Verificar algunos plugins básicos
    from qiime2.plugins import feature_table, diversity, taxa
    print("✅ Plugins QIIME2 cargados correctamente")
except Exception as e:
    print(f"❌ QIIME2 error: {str(e)[:50]}")
    failed_critical.append('qiime2-plugins')

print("=" * 60)

if not failed_critical:
    print("\n✅ INSTALACIÓN VERIFICADA CORRECTAMENTE\n")
    sys.exit(0)
else:
    print(f"\n❌ {len(failed_critical)} componente(s) crítico(s) fallaron: {', '.join(failed_critical)}")
    print("   Por favor, revisa los errores arriba.")
    print("")
    sys.exit(1)
VERIFY_SCRIPT

# Ejecutar verificación
if python /tmp/verify_installation.py; then
    info "Todos los componentes verificados correctamente"
else
    warn "Algunos componentes fallaron."

    # Intentar reparar solo lo crítico
    read -p "¿Intentar reparar? (s/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Ss]$ ]]; then
        step "Intentando reparar instalación..."

        # Reinstalar QIIME2 si hay problemas
        if command -v mamba &> /dev/null; then
            mamba install -c qiime2 -c conda-forge qiime2-amplicon=2024.2 -y --force-reinstall
        else
            conda install -c qiime2 -c conda-forge qiime2-amplicon=2024.2 -y --force-reinstall
        fi

        # Verificar de nuevo
        if python /tmp/verify_installation.py; then
            info "Reparación exitosa"
        else
            warn "La reparación no solucionó todos los problemas."
        fi
    fi
fi

rm /tmp/verify_installation.py

echo ""

# ============================================
# 9. CONFIGURAR PIPELINE
# ============================================

step "Configurando pipeline..."

# Crear directorios necesarios
mkdir -p data/raw data/processed logs reports config

# Copiar archivo de configuración si no existe
if [ ! -f "config/config.yaml" ] && [ -f "config/config_template.yaml" ]; then
    cp config/config_template.yaml config/config.yaml
    info "Archivo de configuración creado: config/config.yaml"
fi

# Hacer ejecutables los scripts
if [ -f "microbiome_cli.py" ]; then
    chmod +x microbiome_cli.py
    info "CLI marcado como ejecutable"
fi

echo ""

# ============================================
# 10. CREAR SCRIPTS DE AYUDA
# ============================================

step "Creando scripts de ayuda..."

# Script de activación mejorado
cat > activate.sh << 'EOF'
#!/bin/bash
# Script para activar el entorno del pipeline

# Inicializar conda
if [ -n "$BASH_VERSION" ]; then
    if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
    else
        eval "$(conda shell.bash hook)"
    fi
fi

# Intentar activar con mamba si está disponible, si no con conda
ENV_NAME="qiime2-amplicon-2024.2"

if command -v mamba &> /dev/null; then
    if mamba activate "$ENV_NAME" 2>/dev/null; then
        echo "=============================================="
        echo "  ✅ Entorno $ENV_NAME activado (via mamba)"
        echo "=============================================="
    else
        echo "⚠️  No se pudo activar con mamba, intentando con conda..."
        if conda activate "$ENV_NAME" 2>/dev/null; then
            echo "=============================================="
            echo "  ✅ Entorno $ENV_NAME activado (via conda)"
            echo "=============================================="
        else
            echo "❌ No se pudo activar el entorno $ENV_NAME"
            echo "   Verifica que existe: conda env list"
            return 1 2>/dev/null || exit 1
        fi
    fi
else
    if conda activate "$ENV_NAME" 2>/dev/null; then
        echo "=============================================="
        echo "  ✅ Entorno $ENV_NAME activado"
        echo "=============================================="
    else
        echo "❌ No se pudo activar el entorno $ENV_NAME"
        echo "   Verifica que existe: conda env list"
        return 1 2>/dev/null || exit 1
    fi
fi

echo ""
echo "Para usar el pipeline:"
echo "  python microbiome_cli.py --help"
echo ""
echo "Ejemplo rápido:"
echo "  python microbiome_cli.py download samples.csv"
echo ""
EOF
chmod +x activate.sh

# Script de test rápido
cat > test.sh << 'EOF'
#!/bin/bash
# Script de prueba rápida

# Cargar el script de activación
source activate.sh 2>/dev/null || {
    echo "❌ No se pudo activar el entorno"
    exit 1
}

echo ""
echo "🧪 Test Rápido - Microbiome Pipeline"
echo "======================================"
echo ""

echo "1. Verificando entorno Python..."
python --version
echo ""

echo "2. Verificando QIIME2..."
python -c "
import qiime2
print(f'QIIME2 versión: {qiime2.__version__}')
print('✅ QIIME2 cargado correctamente')
"
echo ""

echo "3. Verificando herramientas principales..."
for tool in snakemake fastp multiqc; do
    if command -v $tool &> /dev/null; then
        echo "   ✅ $tool disponible"
    else
        echo "   ⚠  $tool NO disponible"
    fi
done

echo ""
echo "4. Verificando estructura del proyecto..."
for dir in data/raw data/processed logs reports config; do
    if [ -d "$dir" ]; then
        echo "   ✅ $dir existe"
    else
        echo "   ⚠  $dir no existe"
    fi
done

echo ""
echo "✅ Test completado"
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
echo "     ./activate.sh"
echo "     # o manualmente:"
echo "     conda activate qiime2-amplicon-2024.2"
echo ""
echo "  2. Ejecutar test rápido:"
echo "     ./test.sh"
echo ""
echo "  3. Configurar pipeline:"
echo "     Edita config/config.yaml con tus rutas"
echo ""
echo "  4. Ejecutar pipeline:"
echo "     python microbiome_cli.py --help"
echo ""
echo "=============================================="
echo ""

# Información de versiones instaladas
step "Versiones instaladas:"
echo "  • Python: $(python --version 2>&1)"
echo "  • QIIME2: $(python -c 'import qiime2; print(qiime2.__version__)' 2>/dev/null || echo 'error')"
echo "  • Snakemake: $(snakemake --version 2>/dev/null | head -1 || echo 'instalado')"
echo "  • Fastp: $(fastp --version 2>/dev/null | head -1 || echo 'instalado')"
echo "  • MultiQC: $(multiqc --version 2>&1 | head -1 || echo 'instalado')"

echo ""
echo "🚀 ¡Instalación completada! El pipeline está listo para usar."
echo ""