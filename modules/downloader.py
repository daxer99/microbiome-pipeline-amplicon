"""
Módulo para descarga de SRA con eliminación automática de archivos SRA
Ahora incluye fallback a ENA cuando SRA Toolkit falla por problemas de SSL
"""
import pandas as pd
import subprocess
from pathlib import Path
import sys
import os
import requests
import time
import shutil
from tqdm import tqdm


def download_sra_from_csv(csv_file, output_dir="data/raw"):
    """Descarga SRA desde archivo CSV, con fallback a ENA si es necesario"""

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

    # Verificar si SRA Toolkit funciona
    sra_works = test_sra_toolkit()

    if not sra_works:
        print("⚠️  SRA Toolkit tiene problemas (posiblemente SSL). Usando ENA...")
        # Descargar desde ENA
        for accession in accessions:
            download_from_ena(accession.strip(), output_dir)
    else:
        # Usar SRA Toolkit como siempre
        for accession in accessions:
            download_single_sra(accession.strip(), output_dir)

    print("✅ Descargas completadas")


def test_sra_toolkit():
    """Prueba si SRA Toolkit funciona correctamente"""
    try:
        # Prueba simple con prefetch --help
        result = subprocess.run(
            ['prefetch', '--version'],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            print("✅ SRA Toolkit funciona correctamente")
            return True
        return False
    except (subprocess.TimeoutExpired, FileNotFoundError, PermissionError) as e:
        print(f"⚠️  SRA Toolkit no funciona: {e}")
        return False
    except Exception as e:
        # Capturar errores de SSL específicamente
        if "SSL" in str(e) or "certificate" in str(e).lower():
            print("⚠️  Problema de SSL con SRA Toolkit")
        return False


def download_single_sra(accession, output_dir):
    """Descarga un solo archivo SRA"""
    print(f"⬇️  Descargando {accession} con SRA Toolkit...")

    try:
        # Crear directorio de salida
        os.makedirs(f"{output_dir}/{accession}", exist_ok=True)

        # Prefetch
        prefetch_cmd = ['prefetch', accession, '-O', output_dir]
        print(f"   Comando: {' '.join(prefetch_cmd)}")

        result = subprocess.run(
            prefetch_cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minutos timeout
        )

        if result.returncode != 0:
            print(f"❌ SRA Toolkit falló para {accession}")
            print(f"   Error: {result.stderr[:200]}")
            print(f"   Intentando con ENA...")
            download_from_ena(accession, output_dir)
            return

        # Verificar si el archivo SRA se descargó
        sra_path = f"{output_dir}/{accession}/{accession}.sra"
        if not os.path.exists(sra_path):
            # Buscar en otros posibles lugares
            possible_paths = [
                f"{output_dir}/{accession}.sra",
                f"{output_dir}/{accession}/{accession}",
                f"{output_dir}/{accession}"
            ]

            for path in possible_paths:
                if os.path.exists(path):
                    sra_path = path
                    break

        if not os.path.exists(sra_path):
            print(f"❌ Archivo SRA no encontrado para {accession}")
            print(f"   Intentando con ENA...")
            download_from_ena(accession, output_dir)
            return

        print(f"✅ Prefetch completado para {accession}")

        # Detectar si es single-end o paired-end
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

    except subprocess.TimeoutExpired:
        print(f"❌ Timeout para {accession}. Intentando con ENA...")
        download_from_ena(accession, output_dir)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error con SRA Toolkit para {accession}: {e}")
        print(f"   Intentando con ENA...")
        download_from_ena(accession, output_dir)
    except Exception as e:
        print(f"❌ Error inesperado para {accession}: {e}")
        print(f"   Intentando con ENA...")
        download_from_ena(accession, output_dir)


def download_from_ena(accession, output_dir):
    """Descarga desde ENA (European Nucleotide Archive)"""
    print(f"🌐 Descargando {accession} desde ENA...")

    # Crear directorio
    sample_dir = Path(output_dir) / accession
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Intentar diferentes métodos de descarga de ENA
    success = False

    # 1: API directa de ENA (más confiable)
    if not success:
        success = download_ena_method1(accession, sample_dir)

    #2: Patrones comunes de URL
    if not success:
        success = download_ena_method2(accession, sample_dir)

    # 3: Google Cloud Archive
    if not success:
        success = download_ena_method3(accession, sample_dir)

    if success:
        print(f"✅ {accession} descargado desde ENA")

        # Crear archivos FASTQ si se descargó un archivo comprimido
        fastq_files = list(sample_dir.glob("*.fastq*"))
        if not fastq_files:
            print(f"⚠️  No se encontraron archivos FASTQ para {accession}")
    else:
        print(f"❌ No se pudo descargar {accession} desde ENA")
        print(f"💡 Puedes intentar manualmente desde:")
        print(f"   https://www.ebi.ac.uk/ena/browser/view/{accession}")


def download_ena_method1(accession, sample_dir):
    """1: Usar API de ENA para obtener URLs de FTP"""
    try:
        # Obtener metadatos de ENA
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport"
        params = {
            'accession': accession,
            'result': 'read_run',
            'fields': 'fastq_ftp,submitted_ftp',
            'format': 'json',
            'download': 'false'
        }

        response = requests.get(url, params=params, timeout=30)

        if response.status_code == 200 and response.text.strip():
            data = response.json()
            if data and isinstance(data, list) and len(data) > 0:
                # Obtener URLs de FTP
                ftp_urls = data[0].get('fastq_ftp', '') or data[0].get('submitted_ftp', '')

                if ftp_urls:
                    urls = ftp_urls.split(';')

                    for i, ftp_url in enumerate(urls):
                        if ftp_url.strip():
                            # Convertir FTP a HTTPS
                            https_url = ftp_url.replace('ftp://', 'https://')
                            filename = f"{accession}_{i + 1}.fastq.gz"
                            filepath = sample_dir / filename

                            if download_file(https_url, filepath):
                                return True
    except Exception as e:
        print(f"   Método 1 falló: {str(e)[:50]}")

    return False


def download_ena_method2(accession, sample_dir):
    """2: Patrones comunes de URL de ENA"""
    try:
        # Patrones de URL comunes
        url_patterns = [
            f"https://www.ebi.ac.uk/ena/browser/api/fastq/{accession}?download=true",
            f"https://www.ebi.ac.uk/ena/data/view/{accession}&display=fastq",
        ]

        # Patrones basados en la estructura típica de ENA
        accession_prefix = accession[:6]
        url_patterns.extend([
            f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession_prefix}/{accession}/{accession}.fastq.gz",
            f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession_prefix}/{accession}/{accession}_1.fastq.gz",
            f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession_prefix}/{accession}/{accession}_2.fastq.gz",
        ])

        for url in url_patterns:
            filename = f"{accession}.fastq.gz"
            if "_1.fastq.gz" in url:
                filename = f"{accession}_1.fastq.gz"
            elif "_2.fastq.gz" in url:
                filename = f"{accession}_2.fastq.gz"

            filepath = sample_dir / filename

            if download_file(url, filepath):
                return True

    except Exception as e:
        print(f"   Método 2 falló: {str(e)[:50]}")

    return False


def download_ena_method3(accession, sample_dir):
    """3: Google Cloud Archive"""
    try:
        url = f"https://storage.googleapis.com/nih-sequence-read-archive/{accession}/{accession}.1"
        filepath = sample_dir / f"{accession}.fastq.gz"

        if download_file(url, filepath):
            return True
    except Exception as e:
        print(f"   Método 3 falló: {str(e)[:50]}")

    return False


def download_file(url, filepath):
    """Descarga un archivo desde una URL"""
    try:
        print(f"   Descargando: {url}")

        # Usar requests con stream para manejar archivos grandes
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()

        # Obtener tamaño del archivo si está disponible
        total_size = int(response.headers.get('content-length', 0))

        # Descargar con barra de progreso
        with open(filepath, 'wb') as f:
            if total_size == 0:
                f.write(response.content)
            else:
                chunk_size = 8192
                with tqdm(total=total_size, unit='B', unit_scale=True,
                          desc=f"   Progreso", ncols=80) as pbar:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))

        # Verificar que el archivo no esté vacío
        if os.path.getsize(filepath) > 1000:
            return True
        else:
            os.remove(filepath)
            return False

    except Exception as e:
        if os.path.exists(filepath):
            os.remove(filepath)
        print(f"   Error descargando: {str(e)[:50]}")
        return False


# Las funciones originales se mantienen igual desde aquí...

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
            subprocess.run([tool, '--version'], capture_output=True, timeout=5)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            missing.append(tool)

    if missing:
        print(f"⚠️  Herramientas SRA faltantes: {', '.join(missing)}")
        print("   Pero no te preocupes, usaremos ENA si es necesario")
        # No retornamos False porque podemos usar ENA como fallback
        return True  # Siempre retornar True para permitir continuar

    return True