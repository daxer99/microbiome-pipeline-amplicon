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
    # Convertir a ruta absoluta para evitar problemas
    output_dir = os.path.abspath(output_dir)

    # Si el directorio de salida existe, eliminarlo completamente
    if os.path.exists(output_dir):
        click.echo(f"🗑️  Eliminando directorio existente: {output_dir}")
        shutil.rmtree(output_dir)

    # NO crear el directorio aquí - PICRUSt2 lo crea automáticamente

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
    # PICRUSt2 guarda los resultados en el subdirectorio 'pathways_out'
    pathways_out_dir = os.path.join(output_dir, 'pathways_out')
    pathway_abundance_tsv_picrust = os.path.join(pathways_out_dir, 'path_abun_unstrat.tsv.gz')

    # Si no existe el archivo .gz, buscar el archivo sin comprimir
    if not os.path.exists(pathway_abundance_tsv_picrust):
        pathway_abundance_tsv_picrust = os.path.join(pathways_out_dir, 'path_abun_unstrat.tsv')

    if not os.path.exists(pathway_abundance_tsv_picrust):
        # Buscar cualquier archivo TSV o TSV.GZ en pathways_out
        if os.path.exists(pathways_out_dir):
            tsv_files = [f for f in os.listdir(pathways_out_dir)
                         if f.endswith('.tsv') or f.endswith('.tsv.gz')]
            if tsv_files:
                pathway_abundance_tsv_picrust = os.path.join(pathways_out_dir, tsv_files[0])
            else:
                all_files = os.listdir(pathways_out_dir)
                raise FileNotFoundError(
                    f"No se encontraron archivos TSV o TSV.GZ en {pathways_out_dir}. "
                    f"Archivos encontrados: {all_files}"
                )
        else:
            all_files = os.listdir(output_dir)
            raise FileNotFoundError(
                f"No se encontró el directorio pathways_out en {output_dir}. "
                f"Archivos encontrados: {all_files}"
            )

    # Convertir el archivo TSV de rutas a formato BIOM y luego a artefacto QIIME2
    pathway_abundance_qza = os.path.join(output_dir, 'pathway_abundance.qza')
    pathway_abundance_tsv_output = os.path.join(output_dir, 'pathway_abundance.tsv')

    try:
        # Leer el TSV de PICRUSt2 (soporta archivos .gz automáticamente)
        click.echo(f"📖 Leyendo archivo: {pathway_abundance_tsv_picrust}")

        # pandas.read_csv detecta automáticamente archivos .gz
        df = pd.read_csv(pathway_abundance_tsv_picrust, sep='\t', index_col=0)
        click.echo(f"✅ Se encontraron {len(df)} rutas metabólicas en {len(df.columns)} muestras")

        # Guardar una copia del TSV en el directorio principal (sin comprimir)
        df.to_csv(pathway_abundance_tsv_output, sep='\t')
        click.echo(f"✅ TSV guardado: {pathway_abundance_tsv_output}")

        # Convertir a formato BIOM
        pathway_abundance_biom_output = os.path.join(output_dir, 'pathway_abundance.biom')
        table = biom.Table(df.values, observation_ids=df.index.tolist(),
                           sample_ids=df.columns.tolist())
        with biom.util.biom_open(pathway_abundance_biom_output, 'w') as f:
            table.to_hdf5(f, "PICRUSt2 pathway abundance")
        click.echo(f"✅ BIOM guardado: {pathway_abundance_biom_output}")

        # Importar a QIIME2
        pathway_abundance = Artifact.import_data('FeatureTable[Frequency]', pathway_abundance_biom_output)
        pathway_abundance.save(pathway_abundance_qza)
        click.echo(f"✅ Artefacto QIIME2 guardado: {pathway_abundance_qza}")

    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        raise RuntimeError(f"Error procesando resultados de PICRUSt2: {str(e)}\n{error_details}")

    return {
        'pathway_abundance_biom': pathway_abundance_biom_output,
        'pathway_abundance_tsv': pathway_abundance_tsv_output,
        'pathway_abundance_qza': pathway_abundance_qza,
        'pathways_out_dir': pathways_out_dir
    }


def normalize_pathway_abundance(pathway_table):
    """Normaliza la abundancia de rutas metabólicas a porcentajes.

    Args:
        pathway_table: Ruta al artefacto QIIME2 de rutas metabólicas.

    Returns:
        DataFrame normalizado.
    """
    if isinstance(pathway_table, str):
        pathway_table = Artifact.load(pathway_table)

    # Obtener la tabla como DataFrame usando view()
    df = pathway_table.view(pd.DataFrame)

    # Normalizar a porcentajes por muestra
    df_normalized = df.div(df.sum(axis=0), axis=1) * 100

    return df_normalized