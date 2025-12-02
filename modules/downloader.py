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
import re


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
    """Descarga un solo archivo SRA (método original)"""
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
    """Descarga desde ENA (European Nucleotide Archive) - Método alternativo"""
    print(f"🌐 Descargando {accession} desde ENA...")

    # Crear directorio
    sample_dir = Path(output_dir) / accession
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Primero, obtener las URLs reales de los archivos FASTQ desde la API de ENA
    fastq_urls = get_ena_fastq_urls(accession)

    if not fastq_urls:
        print(f"❌ No se pudieron obtener URLs FASTQ para {accession} desde ENA")
        return

    print(f"🔗 Encontrados {len(fastq_urls)} archivos FASTQ para {accession}")

    # Descargar cada archivo FASTQ
    success_count = 0
    for i, url in enumerate(fastq_urls):
        # Determinar nombre del archivo
        if len(fastq_urls) == 1:
            filename = f"{accession}.fastq.gz"
        elif len(fastq_urls) == 2:
            if i == 0:
                filename = f"{accession}_1.fastq.gz"
            else:
                filename = f"{accession}_2.fastq.gz"
        else:
            filename = f"{accession}_{i + 1}.fastq.gz"

        filepath = sample_dir / filename

        print(f"   Descargando {filename}...")
        if download_file_with_retry(url, filepath):
            success_count += 1
            print(f"   ✅ {filename} descargado")
        else:
            print(f"   ❌ Error descargando {filename}")

    if success_count > 0:
        print(f"✅ {success_count}/{len(fastq_urls)} archivos descargados para {accession}")

        # Si solo se descargó un archivo pero deberían ser 2, buscar el otro
        if success_count == 1 and len(fastq_urls) == 2:
            print(f"⚠️  Solo se descargó 1 de 2 archivos para {accession}")
            print(f"   Intentando método alternativo...")
            try_download_missing_pair(accession, sample_dir, fastq_urls)
    else:
        print(f"❌ No se pudo descargar ningún archivo para {accession}")
        print(f"💡 Puedes intentar manualmente desde:")
        print(f"   https://www.ebi.ac.uk/ena/browser/view/{accession}")


def get_ena_fastq_urls(accession):
    """Obtiene las URLs reales de los archivos FASTQ desde la API de ENA"""
    urls = []

    try:
        # Método 1: API filereport (más confiable)
        api_url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
        params = {
            'accession': accession,
            'result': 'read_run',
            'fields': 'fastq_ftp,submitted_ftp',
            'format': 'json'
        }

        response = requests.get(api_url, params=params, timeout=30)

        if response.status_code == 200 and response.text.strip():
            try:
                data = response.json()
                if data and isinstance(data, list) and len(data) > 0:
                    # Obtener URLs de FTP
                    ftp_field = data[0].get('fastq_ftp') or data[0].get('submitted_ftp')
                    if ftp_field:
                        ftp_urls = ftp_field.split(';')
                        for ftp_url in ftp_urls:
                            if ftp_url.strip():
                                # Convertir FTP a HTTPS
                                https_url = ftp_url.strip().replace('ftp://', 'https://')
                                urls.append(https_url)
            except ValueError:
                # Si no es JSON, intentar parsear como texto
                pass

        # Si no obtuvo URLs con el método 1, intentar método 2
        if not urls:
            # Método 2: Web scraping de la página de ENA
            page_url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
            response = requests.get(page_url, timeout=30)

            if response.status_code == 200:
                # Buscar URLs de FASTQ en el XML
                content = response.text

                # Patrones para encontrar URLs de FASTQ
                patterns = [
                    r'<FASTQ_FILES>(.*?)</FASTQ_FILES>',
                    r'ftp://ftp\.sra\.ebi\.ac\.uk/vol1/fastq/[^<]*\.fastq\.gz',
                    r'https://[^<]*\.fastq\.gz'
                ]

                for pattern in patterns:
                    matches = re.findall(pattern, content, re.IGNORECASE | re.DOTALL)
                    for match in matches:
                        if 'fastq' in match.lower() and 'ftp.sra.ebi.ac.uk' in match:
                            # Convertir FTP a HTTPS
                            https_url = match.replace('ftp://', 'https://')
                            urls.append(https_url)

        # Método 3: Patrones predefinidos basados en la estructura común de ENA
        if not urls:
            # Estructura común de ENA para paired-end
            accession_6 = accession[:6]
            accession_last3 = accession[-3:] if len(accession) > 6 else accession

            # Patrones para paired-end
            url_patterns = [
                f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession_6}/{accession_last3}/{accession}/{accession}_1.fastq.gz",
                f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession_6}/{accession_last3}/{accession}/{accession}_2.fastq.gz",
                f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession_6}/00{accession_last3[-1]}/{accession}/{accession}_1.fastq.gz",
                f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession_6}/00{accession_last3[-1]}/{accession}/{accession}_2.fastq.gz",
            ]

            # Probar cada patrón
            for pattern in url_patterns:
                # Verificar si la URL existe
                try:
                    head_response = requests.head(pattern, timeout=10, allow_redirects=True)
                    if head_response.status_code == 200:
                        urls.append(pattern)
                except:
                    pass

        # Eliminar duplicados
        urls = list(dict.fromkeys(urls))

        # Si tenemos URLs, devolverlas
        if urls:
            return urls

        # Método 4: Usar el API de descarga directa
        download_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=fastq_ftp&download=true"
        response = requests.get(download_url, timeout=30)

        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            if len(lines) > 1:
                # La segunda línea suele contener las URLs
                fields = lines[1].split('\t')
                if len(fields) > 0:
                    ftp_urls = fields[-1].split(';')
                    for ftp_url in ftp_urls:
                        if ftp_url.strip():
                            https_url = ftp_url.strip().replace('ftp://', 'https://')
                            urls.append(https_url)

    except Exception as e:
        print(f"   Error obteniendo URLs ENA: {str(e)[:100]}")

    return urls


def download_file_with_retry(url, filepath, max_retries=3):
    """Descarga un archivo con reintentos"""
    for attempt in range(max_retries):
        try:
            if download_file(url, filepath):
                return True
            else:
                print(f"   Intento {attempt + 1} falló, reintentando...")
                time.sleep(2)
        except Exception as e:
            print(f"   Error en intento {attempt + 1}: {str(e)[:50]}")
            time.sleep(2)

    return False


def download_file(url, filepath):
    """Descarga un archivo desde una URL"""
    try:
        # Usar requests con stream para manejar archivos grandes
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        }

        response = requests.get(url, stream=True, timeout=60, headers=headers)

        # Si es una redirección, seguirla
        if response.status_code in [301, 302, 303, 307, 308]:
            redirect_url = response.headers.get('Location')
            if redirect_url:
                print(f"   Redireccionando a: {redirect_url[:80]}...")
                response = requests.get(redirect_url, stream=True, timeout=60, headers=headers)

        response.raise_for_status()

        # Verificar que sea un archivo FASTQ (no HTML)
        content_type = response.headers.get('content-type', '').lower()
        if 'html' in content_type or 'text' in content_type:
            # Leer un poco del contenido para verificar
            sample = response.content[:1000].decode('utf-8', errors='ignore')
            if '<html' in sample.lower() or '<!doctype' in sample.lower():
                print(f"   ❌ La URL devuelve HTML, no FASTQ")
                return False

        # Obtener tamaño del archivo si está disponible
        total_size = int(response.headers.get('content-length', 0))

        # Descargar con barra de progreso
        with open(filepath, 'wb') as f:
            if total_size == 0:
                f.write(response.content)
            else:
                chunk_size = 8192
                with tqdm(total=total_size, unit='B', unit_scale=True,
                          desc=f"      Progreso", ncols=80, leave=False) as pbar:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))

        # Verificar que el archivo no esté vacío y sea un archivo FASTQ válido
        file_size = os.path.getsize(filepath)
        if file_size > 1000:
            # Verificar que sea un archivo gzip válido (los FASTQ.gz)
            try:
                import gzip
                with gzip.open(filepath, 'rb') as test_f:
                    test_f.read(100)  # Leer primeros 100 bytes
                return True
            except:
                # Puede que no sea gzip, verificar si es FASTQ plano
                with open(filepath, 'rb') as test_f:
                    header = test_f.read(100)
                    if b'@' in header:  # FASTQ comienza con @
                        return True
                print(f"   ❌ Archivo descargado no es FASTQ válido")
                os.remove(filepath)
                return False
        else:
            os.remove(filepath)
            print(f"   ❌ Archivo demasiado pequeño (posiblemente vacío)")
            return False

    except requests.exceptions.RequestException as e:
        if os.path.exists(filepath):
            os.remove(filepath)
        print(f"   ❌ Error de red: {str(e)[:50]}")
        return False
    except Exception as e:
        if os.path.exists(filepath):
            os.remove(filepath)
        print(f"   ❌ Error descargando: {str(e)[:50]}")
        return False


def try_download_missing_pair(accession, sample_dir, existing_urls):
    """Intenta descargar el archivo FASTQ faltante para paired-end"""
    # Buscar archivos ya descargados
    downloaded_files = list(sample_dir.glob("*.fastq.gz"))

    if len(downloaded_files) == 1:
        downloaded_file = downloaded_files[0].name

        # Determinar cuál falta (R1 o R2)
        if '_1.fastq.gz' in downloaded_file:
            missing_suffix = '_2.fastq.gz'
            existing_suffix = '_1.fastq.gz'
        elif '_2.fastq.gz' in downloaded_file:
            missing_suffix = '_1.fastq.gz'
            existing_suffix = '_2.fastq.gz'
        else:
            # Si no tiene sufijo, asumir que es R1
            missing_suffix = '_2.fastq.gz'
            existing_suffix = ''

        # Intentar construir la URL del archivo faltante
        if existing_urls and len(existing_urls) > 0:
            existing_url = existing_urls[0]
            # Modificar la URL para obtener el otro par
            if '_1.fastq.gz' in existing_url:
                missing_url = existing_url.replace('_1.fastq.gz', '_2.fastq.gz')
            elif '_2.fastq.gz' in existing_url:
                missing_url = existing_url.replace('_2.fastq.gz', '_1.fastq.gz')
            else:
                # Si la URL no tiene sufijo, intentar patrón común
                accession_6 = accession[:6]
                accession_last3 = accession[-3:] if len(accession) > 6 else accession

                if missing_suffix == '_2.fastq.gz':
                    missing_url = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession_6}/{accession_last3}/{accession}/{accession}_2.fastq.gz"
                else:
                    missing_url = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{accession_6}/{accession_last3}/{accession}/{accession}_1.fastq.gz"

            missing_filename = f"{accession}{missing_suffix}"
            missing_filepath = sample_dir / missing_filename

            print(f"   Intentando descargar {missing_filename}...")
            if download_file_with_retry(missing_url, missing_filepath):
                print(f"   ✅ {missing_filename} descargado")
            else:
                print(f"   ❌ No se pudo descargar {missing_filename}")


# Las funciones restantes se mantienen igual...

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