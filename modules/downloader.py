"""
Módulo para descarga de SRA con eliminación automática de archivos SRA
"""
import pandas as pd
import subprocess
from pathlib import Path
import sys
import os


def download_sra_from_csv(csv_file, output_dir="data/raw"):
    """Descarga SRA desde archivo CSV"""

    # Leer CSV
    df = pd.read_csv(csv_file)

    # Buscar columna de accessions
    accession_col = None
    for col in df.columns:
        if any(x in col.lower() for x in ['accession', 'sra', 'run']):
            accession_col = col
            break

    if not accession_col:
        print("❌ No se encontró columna de accessions")
        sys.exit(1)

    accessions = df[accession_col].dropna().unique()
    print(f"📥 Descargando {len(accessions)} muestras...")

    # Descargar cada accession
    for accession in accessions:
        download_single_sra(accession.strip(), output_dir)

    print("✅ Descargas completadas")


def download_single_sra(accession, output_dir):
    """Descarga un solo archivo SRA y detecta si es SE o PE"""
    print(f"⬇️  Descargando {accession}...")

    try:
        # Prefetch
        subprocess.run([
            'prefetch', accession, '-O', output_dir
        ], check=True, capture_output=True)

        # Detectar si es single-end o paired-end
        sra_path = f"{output_dir}/{accession}/{accession}.sra"
        is_paired = detect_paired_end(sra_path)

        print(f"🔍 {accession} detectado como {'Paired-End' if is_paired else 'Single-End'}")

        # Convertir a FASTQ según el tipo
        if is_paired:
            convert_paired_end(sra_path, output_dir, accession)
        else:
            convert_single_end(sra_path, output_dir, accession)

        # Eliminar archivos SRA después de la conversión
        cleanup_sra_files(output_dir, accession)

        print(f"✅ {accession} descargado, convertido y limpiado")

    except subprocess.CalledProcessError as e:
        print(f"❌ Error con {accession}: {e}")


def detect_paired_end(sra_file):
    """Detecta si el archivo SRA es paired-end"""
    try:
        # Usar sra-stat para ver la estructura del archivo
        result = subprocess.run([
            'sra-stat', '--quick', '--xml', sra_file
        ], capture_output=True, text=True, check=True)

        # Buscar indicios de paired-end en el XML
        if 'read="2"' in result.stdout or 'R2' in result.stdout:
            return True

        # Si no encuentra evidencia de paired-end, asumir single-end
        return False

    except subprocess.CalledProcessError:
        # Si sra-stat falla, usar método alternativo
        return detect_paired_end_alt(sra_file)


def detect_paired_end_alt(sra_file):
    """Método alternativo para detectar paired-end"""
    try:
        # Usar fasterq-dump en modo dry-run para ver la estructura
        result = subprocess.run([
            'fasterq-dump', '--split-3', '--dry-run', sra_file
        ], capture_output=True, text=True)

        # Si muestra dos archivos (_1.fastq y _2.fastq), es paired-end
        output = result.stdout + result.stderr
        if '_1.fastq' in output and '_2.fastq' in output:
            return True

        return False

    except Exception:
        # Por defecto, asumir single-end
        print(f"⚠️  No se pudo determinar el tipo de {sra_file}, asumiendo Single-End")
        return False


def convert_single_end(sra_file, output_dir, accession):
    """Convierte SRA single-end a FASTQ"""
    subprocess.run([
        'fasterq-dump',
        sra_file,
        '--outdir', f"{output_dir}/{accession}",
        '--skip-technical',
        '--threads', '2'
    ], check=True)


def convert_paired_end(sra_file, output_dir, accession):
    """Convierte SRA paired-end a FASTQ"""
    subprocess.run([
        'fasterq-dump',
        sra_file,
        '--outdir', f"{output_dir}/{accession}",
        '--split-files',
        '--threads', '2'
    ], check=True)


def cleanup_sra_files(output_dir, accession):
    """Elimina archivos SRA después de la conversión a FASTQ"""
    accession_dir = Path(output_dir) / accession

    # Buscar y eliminar archivos SRA
    sra_files = list(accession_dir.glob("*.sra"))
    for sra_file in sra_files:
        try:
            sra_file.unlink()
            print(f"🗑️  Eliminado: {sra_file}")
        except Exception as e:
            print(f"⚠️  No se pudo eliminar {sra_file}: {e}")

    # También eliminar archivos intermedios como .csi si existen
    temp_files = list(accession_dir.glob("*.csi")) + list(accession_dir.glob("*.vdbcache"))
    for temp_file in temp_files:
        try:
            temp_file.unlink()
        except Exception:
            pass  # Ignorar errores en archivos temporales


def check_dependencies():
    """Verifica que las herramientas necesarias estén instaladas"""
    tools = ['prefetch', 'fasterq-dump', 'sra-stat']
    missing = []

    for tool in tools:
        try:
            subprocess.run([tool, '--version'], capture_output=True)
        except FileNotFoundError:
            missing.append(tool)

    if missing:
        print(f"❌ Herramientas faltantes: {', '.join(missing)}")
        print("Instala SRA Toolkit: https://github.com/ncbi/sra-tools")
        return False

    return True


def get_download_stats(output_dir):
    """Muestra estadísticas de los archivos descargados"""
    output_path = Path(output_dir)

    if not output_path.exists():
        print("❌ El directorio de salida no existe")
        return

    fastq_files = list(output_path.rglob("*.fastq"))
    sra_files = list(output_path.rglob("*.sra"))

    print(f"\n📊 Estadísticas de descarga:")
    print(f"   • Archivos FASTQ: {len(fastq_files)}")
    print(f"   • Archivos SRA restantes: {len(sra_files)}")

    # Tamaño total de FASTQ
    total_size = sum(f.stat().st_size for f in fastq_files)
    print(f"   • Espacio usado por FASTQ: {total_size / (1024 ** 3):.2f} GB")