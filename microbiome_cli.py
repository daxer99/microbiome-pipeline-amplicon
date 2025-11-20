#!/usr/bin/env python3
"""
CLI Principal - Análisis de Microbioma 16S

Pipeline completo para análisis de datos de microbioma 16S:
- Descarga de secuencias SRA desde CSV
- Creación de manifiestos para QIIME2
- Importación a QIIME2
- Control de calidad y filtrado
- Denoising con Deblur

Ejemplos de uso:
  python microbiome_cli.py download samples.csv
  python microbiome_cli.py quality-control demux.qza
  python microbiome_cli.py run-deblur demux.qza --trim-length 250
"""
import click
from modules.downloader import download_sra_from_csv, check_dependencies
from modules.qiime2_utils import create_fasta_manifest, import_to_qiime2, check_qiime2_installation
from modules.quality_control import QualityControl
from modules.denoiser import Denoiser


@click.group(invoke_without_command=True)
@click.pass_context
@click.option('--version', '-v', is_flag=True, help='Mostrar versión')
def cli(ctx, version):
    """Herramienta de análisis de microbioma 16S

    Un pipeline completo para el análisis de datos de amplicón 16S
    desde la descarga de secuencias hasta la obtención de ASVs.
    """
    if version:
        click.echo("Microbiome Pipeline 16S v1.0.0")
        return

    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        click.echo("\n📋 Comandos disponibles:")
        click.echo("  download         Descargar secuencias SRA desde CSV")
        click.echo("  create-manifest  Crear archivo de manifiesto para QIIME2")
        click.echo("  import-qiime2    Importar datos a QIIME2")
        click.echo("  quality-control  Control de calidad completo")
        click.echo("  run-deblur       Denoising con Deblur")
        click.echo("\n💡 Usa 'microbiome_cli.py COMANDO --help' para ayuda específica")


@cli.command()
@click.argument('csv_file', type=click.Path(exists=True))
@click.option('--output-dir', default='data/raw', help='Directorio de salida para los FASTQ')
@click.option('--accession-col', default='run_accession',
              help='Nombre de la columna con los accessions (por defecto: run_accession)')
def download(csv_file, output_dir, accession_col):
    """Descargar archivos SRA desde un archivo CSV

    CSV_FILE: Ruta al archivo CSV que contiene los accessions SRA

    Ejemplos:
      microbiome_cli.py download samples.csv
      microbiome_cli.py download samples.csv --output-dir my_data
      microbiome_cli.py download samples.csv --accession-col sample_id
    """
    if not check_dependencies():
        return

    click.echo(f"📥 Descargando secuencias desde: {csv_file}")
    click.echo(f"📁 Directorio de salida: {output_dir}")
    click.echo(f"🔤 Columna de accessions: {accession_col}")

    download_sra_from_csv(csv_file, output_dir)


@cli.command()
@click.option('--input-dir', default='data/raw',
              help='Directorio con archivos FASTQ (por defecto: data/raw)')
@click.option('--output-file', default='fasta_manifest.csv',
              help='Archivo de manifiesto de salida (por defecto: fasta_manifest.csv)')
def create_manifest(input_dir, output_file):
    """Crear archivo de manifiesto para importación en QIIME2

    Busca automáticamente archivos FASTQ en el directorio de entrada
    y genera un manifiesto en formato CSV compatible con QIIME2.

    Ejemplos:
      microbiome_cli.py create-manifest
      microbiome_cli.py create-manifest --input-dir my_data
      microbiome_cli.py create-manifest --output-file my_manifest.csv
    """
    click.echo(f"🔍 Buscando FASTQ en: {input_dir}")
    click.echo(f"📄 Creando manifiesto: {output_file}")

    create_fasta_manifest(input_dir, output_file)


@cli.command()
@click.argument('manifest_file', type=click.Path(exists=True))
@click.option('--output-dir', default='data/qiime2',
              help='Directorio de salida para artefactos QIIME2 (por defecto: data/qiime2)')
def import_qiime2(manifest_file, output_dir):
    """Importar datos a QIIME2 desde archivo de manifiesto

    MANIFEST_FILE: Ruta al archivo de manifiesto CSV

    Convierte los archivos FASTQ en un artefacto QIIME2 (.qza)
    para análisis posteriores.

    Ejemplos:
      microbiome_cli.py import-qiime2 manifest.csv
      microbiome_cli.py import-qiime2 manifest.csv --output-dir qiime_data
    """
    if not check_qiime2_installation():
        return

    click.echo(f"📦 Importando desde: {manifest_file}")
    click.echo(f"📁 Directorio de salida: {output_dir}")

    import_to_qiime2(manifest_file, output_dir)


@cli.command()
@click.argument('demux_file', type=click.Path(exists=True))
@click.option('--output-dir', default='results/quality_control',
              help='Directorio de salida (por defecto: results/quality_control)')
@click.option('--min-quality', default=20,
              help='Calidad mínima para filtrado (por defecto: 20)')
def quality_control(demux_file, output_dir, min_quality):
    """Control de calidad completo: reportes, gráficos y filtrado

    DEMUX_FILE: Ruta al artefacto QIIME2 con secuencias demultiplexadas (.qza)

    Genera:
      - Reporte de calidad QIIME2 (.qzv)
      - Gráficos de perfil de calidad (.png)
      - Secuencias filtradas (.qza)

    Ejemplos:
      microbiome_cli.py quality-control demux.qza
      microbiome_cli.py quality-control demux.qza --min-quality 25
      microbiome_cli.py quality-control demux.qza --output-dir my_qc
    """
    click.echo(f"🎯 Analizando: {demux_file}")
    click.echo(f"📁 Directorio de salida: {output_dir}")
    click.echo(f"📊 Calidad mínima: {min_quality}")

    qc = QualityControl(demux_file)
    results = qc.run_quality_control(output_dir, min_quality)

    if results and results['filtered_seqs']:
        click.echo(f"✅ Secuencias filtradas: {results['filtered_seqs']}")


@cli.command()
@click.argument('demux_file', type=click.Path(exists=True))
@click.option('--output-dir', default='results/deblur',
              help='Directorio de salida (por defecto: results/deblur)')
@click.option('--left-trim-len', default=0,
              help='Longitud de trim de inicio (por defecto: 0)')
@click.option('--trim-length', default=250,
              help='Longitud de trim final (por defecto: 250)')
@click.option('--min-reads', default=10,
              help='Mínimo de lecturas por muestra (por defecto: 10)')
@click.option('--min-size', default=2,
              help='Mínimo de tamaño para filtrado (por defecto: 2)')
@click.option('--jobs-to-start', default=8,
              help='Número de CPUs para procesar (por defecto: 8)')
def run_deblur(demux_file, output_dir, left_trim_len, trim_length, min_reads, min_size, jobs_to_start):
    """Ejecutar Deblur para denoising y obtención de ASVs

    DEMUX_FILE: Ruta al artefacto QIIME2 con secuencias demultiplexadas (.qza)

    Deblur es un método rápido para obtener ASVs (Amplicon Sequence Variants)
    mediante corrección de errores.

    Ejemplos:
      microbiome_cli.py run-deblur demux.qza
      microbiome_cli.py run-deblur demux.qza --trim-length 200
      microbiome_cli.py run-deblur demux.qza --jobs-to-start 4
    """
    click.echo(f"🧹 Ejecutando Deblur en: {demux_file}")
    click.echo(f"📁 Directorio de salida: {output_dir}")
    click.echo(f"✂️  Trim inicial: {left_trim_len}")
    click.echo(f"✂️  Trim final: {trim_length}")
    click.echo(f"📊 Mínimo de lecturas: {min_reads}")
    click.echo(f"🔢 Mínimo de tamaño: {min_size}")
    click.echo(f"⚡ CPUs: {jobs_to_start}")

    denoiser = Denoiser(demux_file)
    result = denoiser.run_deblur(
        output_dir=output_dir,
        left_trim_len=left_trim_len,
        trim_length=trim_length,
        min_reads=min_reads,
        min_size=min_size,
        jobs_to_start=jobs_to_start
    )

    if result:
        click.echo("🎉 Proceso de denoising con Deblur completado exitosamente!")


if __name__ == '__main__':
    cli()