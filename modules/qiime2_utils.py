"""
Utilidades para QIIME2 - Generación de manifest files
"""
import pandas as pd
from pathlib import Path
import os
import subprocess
from qiime2 import Artifact


def create_fasta_manifest(input_dir, output_file="fasta_manifest.csv"):
    """
    Crea un archivo de manifiesto para QIIME2 a partir de archivos FASTQ
    """
    input_path = Path(input_dir)
    manifest_data = []

    print(f"🔍 Buscando archivos FASTQ en {input_dir}...")

    # Buscar todas las carpetas de muestras
    sample_dirs = [d for d in input_path.iterdir() if d.is_dir()]

    if not sample_dirs:
        print("❌ No se encontraron carpetas de muestras")
        return None

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name

        # Buscar archivos FASTQ en la carpeta de la muestra
        fastq_files = list(sample_dir.glob("*.fastq")) + list(sample_dir.glob("*.fastq.gz"))

        if not fastq_files:
            print(f"⚠️  No se encontraron archivos FASTQ en {sample_dir}")
            continue

        # Determinar si es single-end o paired-end
        if len(fastq_files) == 1:
            # Single-end
            manifest_data.append({
                'sample-id': sample_id,
                'absolute-filepath': str(fastq_files[0].resolve()),
                'direction': 'forward'
            })
            print(f"📄 {sample_id}: Single-end -> {fastq_files[0].name}")

        elif len(fastq_files) == 2:
            # Paired-end - identificar forward y reverse
            forward, reverse = identify_reads(fastq_files, sample_id)

            if forward and reverse:
                manifest_data.append({
                    'sample-id': sample_id,
                    'absolute-filepath': str(forward.resolve()),
                    'direction': 'forward'
                })
                manifest_data.append({
                    'sample-id': sample_id,
                    'absolute-filepath': str(reverse.resolve()),
                    'direction': 'reverse'
                })
                print(f"📄 {sample_id}: Paired-end -> {forward.name}, {reverse.name}")
            else:
                print(f"⚠️  No se pudieron identificar reads forward/reverse para {sample_id}")
        else:
            print(f"⚠️  Número inesperado de archivos FASTQ en {sample_dir}: {len(fastq_files)}")

    if not manifest_data:
        print("❌ No se encontraron datos válidos para el manifiesto")
        return None

    # Crear DataFrame y guardar en formato CSV
    df = pd.DataFrame(manifest_data)
    df.to_csv(output_file, index=False)

    print(f"✅ Manifest file creado: {output_file}")
    print(f"   • Muestras procesadas: {len(set([d['sample-id'] for d in manifest_data]))}")
    print(f"   • Entradas en el manifest: {len(manifest_data)}")

    return output_file


def identify_reads(fastq_files, sample_id):
    """
    Identifica qué archivo es forward y cuál es reverse
    """
    # Patrones comunes para identificar forward/reverse
    forward_patterns = ['_1', '_R1', '_R1_', '.1.', '_forward', '_F', 'R1']
    reverse_patterns = ['_2', '_R2', '_R2_', '.2.', '_reverse', '_R', 'R2']

    forward_candidates = []
    reverse_candidates = []

    for fq_file in fastq_files:
        filename = fq_file.name

        # Verificar patrones forward
        if any(pattern in filename for pattern in forward_patterns):
            forward_candidates.append(fq_file)
        # Verificar patrones reverse
        elif any(pattern in filename for pattern in reverse_patterns):
            reverse_candidates.append(fq_file)

    # Si encontramos exactamente uno de cada, perfecto
    if len(forward_candidates) == 1 and len(reverse_candidates) == 1:
        return forward_candidates[0], reverse_candidates[0]

    # Si no, intentar por orden alfabético
    sorted_files = sorted(fastq_files)
    if len(sorted_files) == 2:
        print(
            f"ℹ️  Usando orden alfabético para {sample_id}: {sorted_files[0].name} -> forward, {sorted_files[1].name} -> reverse")
        return sorted_files[0], sorted_files[1]

    print(f"❌ No se pudieron identificar reads para {sample_id}")
    return None, None


def import_to_qiime2(manifest_file, output_dir="data/qiime2"):
    """
    Importa los archivos FASTQ a un artefacto QIIME2 usando la API directa
    """
    manifest_path = Path(manifest_file)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if not manifest_path.exists():
        print(f"❌ El archivo de manifiesto no existe: {manifest_file}")
        return None

    # Determinar si es single-end o paired-end
    df = pd.read_csv(manifest_file)
    directions = df['direction'].unique()

    print(f"📊 Archivo de manifest detectado con {len(df)} entradas")
    print(f"   • Direcciones encontradas: {list(directions)}")

    try:
        if 'reverse' in directions:
            # Paired-end
            print("🔍 Tipo detectado: Paired-End")
            paired_end_demux = Artifact.import_data(
                'SampleData[PairedEndSequencesWithQuality]',
                manifest_path,
                view_type='PairedEndFastqManifestPhred33'
            )
            output_file = output_path / "paired_end_demux.qza"
            paired_end_demux.save(str(output_file))
            print(f"✅ Artefacto QIIME2 creado: {output_file}")
            return str(output_file)
        else:
            # Single-end
            print("🔍 Tipo detectado: Single-End")
            single_end_demux = Artifact.import_data(
                'SampleData[SequencesWithQuality]',
                manifest_path,
                view_type='SingleEndFastqManifestPhred33V2'
            )
            output_file = output_path / "single_end_demux.qza"
            single_end_demux.save(str(output_file))
            print(f"✅ Artefacto QIIME2 creado: {output_file}")
            return str(output_file)

    except Exception as e:
        print(f"❌ Error importando con la API de QIIME2: {e}")
        return None


def check_qiime2_installation():
    """Verifica que QIIME2 esté instalado y disponible"""
    try:
        # Solo verificamos que podemos importar Artifact
        from qiime2 import Artifact
        print("✅ API de QIIME2 detectada")
        return True
    except ImportError as e:
        print(f"❌ No se puede importar QIIME2: {e}")
        print("   Asegúrate de estar en el entorno de QIIME2 correcto")
        return False