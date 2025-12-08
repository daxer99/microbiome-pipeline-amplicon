#!/bin/bash
# install_from_environment.sh
# Instalación usando el environment.yml existente

echo "=============================================="
echo "   Instalación desde environment.yml"
echo "=============================================="

# 1. Inicializar Miniforge3
source ~/miniforge3/etc/profile.d/conda.sh

# 2. Instalar mamba si no está
if ! command -v mamba &> /dev/null; then
    echo "Instalando mamba..."
    conda install -n base -c conda-forge mamba -y
fi

# 3. Limpiar entorno anterior
echo "Limpiando entorno anterior..."
conda deactivate 2>/dev/null || true
conda env remove -n qiime2-amplicon-2024.2 -y 2>/dev/null || true
rm -rf ~/miniforge3/envs/qiime2-amplicon-2024.2 2>/dev/null || true

# 4. Verificar que environment.yml existe
if [ ! -f "environment.yml" ]; then
    echo "❌ Error: No se encuentra environment.yml"
    exit 1
fi

echo "Usando environment.yml para crear el entorno..."
echo "Esto puede tardar 20-30 minutos..."

# 5. Crear entorno desde environment.yml
mamba env create -n qiime2-amplicon-2024.2 -f environment.yml

if [ $? -eq 0 ]; then
    echo "✅ Entorno creado exitosamente"
else
    echo "❌ Error al crear entorno. Intentando con conda..."
    conda env create -n qiime2-amplicon-2024.2 -f environment.yml || {
        echo "❌ Error crítico: No se pudo crear el entorno"
        exit 1
    }
fi

# 6. Activar entorno
echo "Activando entorno..."
mamba activate qiime2-amplicon-2024.2 || conda activate qiime2-amplicon-2024.2

# 7. Instalar dependencias adicionales del pipeline
echo "Instalando dependencias adicionales..."

# Snakemake y herramientas específicas del pipeline
mamba install -c conda-forge -c bioconda \
    snakemake-minimal=7.32.4 \
    multiqc=1.17 \
    fastp=0.23.4 \
    sra-tools \
    -y

# 8. Instalar dokdo
echo "Instalando dokdo..."
pip install git+https://github.com/sbslee/dokdo.git 2>/dev/null || {
    echo "Instalando dokdo con --no-deps..."
    pip install --no-deps git+https://github.com/sbslee/dokdo.git
}

# 9. Verificar instalación
echo "Verificando instalación..."
echo ""
echo "=== VERSIONES INSTALADAS ==="
python -c "import qiime2; print(f'QIIME2: {qiime2.__version__}')" 2>/dev/null || echo "QIIME2: ❌ No disponible"
snakemake --version 2>/dev/null | head -1 || echo "Snakemake: ❌ No disponible"
fastp --version 2>/dev/null | head -1 || echo "Fastp: ❌ No disponible"
multiqc --version 2>&1 | head -1 || echo "MultiQC: ❌ No disponible"

# 10. Configurar directorios
echo ""
echo "Configurando estructura del pipeline..."
mkdir -p data/raw data/processed logs reports config scripts

if [ -f "config/config_template.yaml" ] && [ ! -f "config/config.yaml" ]; then
    cp config/config_template.yaml config/config.yaml
    echo "✅ Archivo de configuración creado: config/config.yaml"
fi

# 11. Crear scripts de ayuda
echo "Creando scripts de ayuda..."

cat > activate_microbiome.sh << 'EOF'
#!/bin/bash
# Script para activar el entorno del pipeline

# Inicializar conda
source ~/miniforge3/etc/profile.d/conda.sh

# Activar entorno
mamba activate qiime2-amplicon-2024.2 2>/dev/null || conda activate qiime2-amplicon-2024.2

if [[ "$CONDA_DEFAULT_ENV" == "qiime2-amplicon-2024.2" ]]; then
    echo "=============================================="
    echo "  ✅ Entorno qiime2-amplicon-2024.2 activado"
    echo "=============================================="
    echo ""
    echo "Herramientas disponibles:"
    echo "  • qiime --help"
    echo "  • snakemake --help"
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

cat > test_installation.sh << 'EOF'
#!/bin/bash
echo "🧪 Probando instalación del Microbiome Pipeline"
echo "=============================================="

# Activar entorno
source activate_microbiome.sh 2>/dev/null || {
    echo "❌ No se pudo activar el entorno"
    exit 1
}

echo ""
echo "1. Verificando herramientas principales:"
tools=("qiime" "snakemake" "fastp" "multiqc" "vsearch" "cutadapt")
for tool in "${tools[@]}"; do
    if command -v $tool &> /dev/null; then
        echo "   ✅ $tool disponible"
    else
        echo "   ❌ $tool NO disponible"
    fi
done

echo ""
echo "2. Verificando Python y módulos:"
python -c "
modules = ['qiime2', 'pandas', 'numpy', 'biom']
for mod in modules:
    try:
        __import__(mod)
        print(f'   ✅ {mod}')
    except ImportError:
        print(f'   ❌ {mod}')
"

echo ""
echo "✅ Prueba completada"
EOF

chmod +x test_installation.sh

# 12. Final
echo ""
echo "=============================================="
echo "   ✅ INSTALACIÓN COMPLETADA"
echo "=============================================="
echo ""
echo "Para usar el pipeline:"
echo ""
echo "  1. Activar el entorno:"
echo "     source activate_microbiome.sh"
echo ""
echo "  2. Probar la instalación:"
echo "     ./test_installation.sh"
echo ""
echo "  3. Configurar el pipeline:"
echo "     Edita config/config.yaml con tus parámetros"
echo ""
echo "  4. Ejecutar el pipeline:"
echo "     python microbiome_cli.py --help"
echo ""