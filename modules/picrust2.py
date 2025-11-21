# modules/picrust2.py
import os
import subprocess
import tempfile
import pathlib
import shutil
from qiime2 import Artifact
import pandas as pd
import biom
from qiime2.plugins.feature_table.methods import filter_features
import click


def check_picrust2_installation():
    """Verifica si PICRUSt2 está instalado correctamente."""
    try:
        result = subprocess.run(['picrust2_pipeline.py', '--version'],
                                capture_output=True, text=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def run_picrust2(table, rep_seqs, output_dir, threads=1):
    """Ejecuta PICRUSt2 para inferir rutas metabólicas.

    Args:
        table: Ruta al artefacto QIIME2 de la tabla de características (.qza).
        rep_seqs: Ruta al artefacto QIIME2 de secuencias representativas (.qza).
        output_dir: Directorio de salida para los resultados.
        threads: Número de hilos a usar.

    Returns:
        dict: Rutas a los archivos de salida.
    """
    # Si el directorio de salida existe, eliminarlo
    if os.path.exists(output_dir):
        click.echo(f"🗑️  Eliminando directorio existente: {output_dir}")
        shutil.rmtree(output_dir)

    os.makedirs(output_dir, exist_ok=True)

    # Verificar que PICRUSt2 esté instalado
    if not check_picrust2_installation():
        raise RuntimeError(
            "PICRUSt2 no está instalado o no está en el PATH. "
            "Instálalo con: conda install -c bioconda picrust2"
        )

    # Cargar artefactos
    table_artifact = Artifact.load(table)
    rep_seqs_artifact = Artifact.load(rep_seqs)

    # Exportar a archivos temporales
    with tempfile.TemporaryDirectory() as tmpdir:
        # Exportar tabla a formato BIOM
        table_export_dir = os.path.join(tmpdir, 'table_export')
        table_artifact.export_data(table_export_dir)

        # Encontrar el archivo BIOM (puede tener diferentes nombres)
        biom_files = list(pathlib.Path(table_export_dir).glob('*.biom'))
        if not biom_files:
            raise FileNotFoundError(f"No se encontró archivo BIOM en {table_export_dir}")
        table_biom_path = str(biom_files[0])

        # Exportar secuencias a FASTA
        seqs_export_dir = os.path.join(tmpdir, 'seqs_export')
        rep_seqs_artifact.export_data(seqs_export_dir)

        # Encontrar el archivo FASTA (puede tener diferentes nombres)
        fasta_files = list(pathlib.Path(seqs_export_dir).glob('*.fasta'))
        if not fasta_files:
            # Intentar con extensión .fna
            fasta_files = list(pathlib.Path(seqs_export_dir).glob('*.fna'))
        if not fasta_files:
            raise FileNotFoundError(f"No se encontró archivo FASTA en {seqs_export_dir}")
        rep_seqs_fasta_path = str(fasta_files[0])

        # Verificar que los archivos se crearon correctamente
        if not os.path.exists(table_biom_path):
            raise FileNotFoundError(f"No se pudo exportar la tabla a: {table_biom_path}")
        if not os.path.exists(rep_seqs_fasta_path):
            raise FileNotFoundError(f"No se pudo exportar las secuencias a: {rep_seqs_fasta_path}")

        click.echo(f"📁 Archivo BIOM: {table_biom_path}")
        click.echo(f"📁 Archivo FASTA: {rep_seqs_fasta_path}")

        # Ejecutar PICRUSt2
        cmd = [
            'picrust2_pipeline.py',
            '-s', rep_seqs_fasta_path,
            '-i', table_biom_path,
            '-o', output_dir,
            '--processes', str(threads),
            '--verbose'  # Añadir verbose para más información
        ]

        try:
            click.echo("🚀 Ejecutando PICRUSt2...")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            if result.stdout:
                click.echo(f"📝 PICRUSt2 stdout: {result.stdout}")
        except subprocess.CalledProcessError as e:
            error_msg = f"PICRUSt2 failed with exit code {e.returncode}:\n"
            if e.stdout:
                error_msg += f"STDOUT: {e.stdout}\n"
            if e.stderr:
                error_msg += f"STDERR: {e.stderr}"
            raise RuntimeError(error_msg)

    # Verificar que los archivos de salida se crearon
    pathway_abundance_biom = os.path.join(output_dir, 'pathway_abundance.biom')
    if not os.path.exists(pathway_abundance_biom):
        # Buscar otros posibles nombres de archivo
        biom_files = [f for f in os.listdir(output_dir) if f.endswith('.biom')]
        if biom_files:
            pathway_abundance_biom = os.path.join(output_dir, biom_files[0])
        else:
            # Listar todos los archivos en el directorio para diagnóstico
            all_files = os.listdir(output_dir)
            raise FileNotFoundError(
                f"No se encontró el archivo pathway_abundance.biom en {output_dir}. "
                f"Archivos encontrados: {all_files}"
            )

    # Convertir el archivo BIOM de rutas a artefacto QIIME2
    pathway_abundance_qza = os.path.join(output_dir, 'pathway_abundance.qza')
    try:
        pathway_abundance = Artifact.import_data('FeatureTable[Frequency]', pathway_abundance_biom)
        pathway_abundance.save(pathway_abundance_qza)
    except Exception as e:
        raise RuntimeError(f"Error importing BIOM to QIIME2: {str(e)}")

    # También generar el archivo TSV
    pathway_abundance_tsv = os.path.join(output_dir, 'pathway_abundance.tsv')
    try:
        table_biom = biom.load_table(pathway_abundance_biom)
        with open(pathway_abundance_tsv, 'w') as f:
            table_biom.to_tsv(header_key='KEGG_Pathways', header_value='KEGG_Pathways', direct_io=f)
    except Exception as e:
        print(f"Warning: Could not create TSV file: {str(e)}")

    return {
        'pathway_abundance_biom': pathway_abundance_biom,
        'pathway_abundance_tsv': pathway_abundance_tsv,
        'pathway_abundance_qza': pathway_abundance_qza
    }


def filter_low_abundance_pathways(pathway_table, min_abundance=0.001):
    """Filtra rutas metabólicas de baja abundancia.

    Args:
        pathway_table: Ruta al artefacto QIIME2 de rutas metabólicas.
        min_abundance: Abundancia mínima para mantener una ruta.

    Returns:
        Artefacto QIIME2 filtrado.
    """
    if isinstance(pathway_table, str):
        pathway_table = Artifact.load(pathway_table)

    # Filtrar características con baja abundancia
    filtered_table = filter_features(
        table=pathway_table,
        min_frequency=min_abundance
    )

    return filtered_table.filtered_table


def normalize_pathway_abundance(pathway_table):
    """Normaliza la abundancia de rutas metabólicas a porcentajes.

    Args:
        pathway_table: Ruta al artefacto QIIME2 de rutas metabólicas.

    Returns:
        DataFrame normalizado.
    """
    if isinstance(pathway_table, str):
        pathway_table = Artifact.load(pathway_table)

    # Exportar a DataFrame temporal
    with tempfile.TemporaryDirectory() as tmpdir:
        pathway_table.export_data(tmpdir)
        data_dir_fp = pathlib.Path(tmpdir)
        csv_files = list(data_dir_fp.glob('*.tsv'))
        if not csv_files:
            csv_files = list(data_dir_fp.glob('*.csv'))
        if not csv_files:
            raise FileNotFoundError(f"No se encontraron archivos TSV o CSV en {tmpdir}")

        csv_file = csv_files[0]
        df = pd.read_csv(csv_file, sep='\t', index_col=0)

    # Normalizar a porcentajes por muestra
    df_normalized = df.div(df.sum(axis=0), axis=1) * 100

    return df_normalized