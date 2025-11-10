#!/usr/bin/env python3
"""
CLI Principal - Análisis de Microbioma 16S
"""
import click
from modules.downloader import download_sra_from_csv, check_dependencies
from modules.qiime2_utils import create_fasta_manifest, import_to_qiime2, check_qiime2_installation
from modules.quality_control import QualityControl
from modules.denoiser import Denoiser


@click.group()
def cli():
    """Herramienta de análisis de microbioma 16S"""
    pass


@cli.command()
@click.argument('csv_file')
@click.option('--output-dir', default='data/raw', help='Directorio de salida')
def download(csv_file, output_dir):
    """Descargar SRA desde CSV"""
    if not check_dependencies():
        return
    download_sra_from_csv(csv_file, output_dir)


@cli.command()
@click.option('--input-dir', default='data/raw', help='Directorio con FASTQ')
@click.option('--output-file', default='fasta_manifest.csv', help='Archivo de manifiesto de salida')
def create_manifest(input_dir, output_file):
    """Crear archivo de manifiesto para QIIME2"""
    create_fasta_manifest(input_dir, output_file)


@cli.command()
@click.argument('manifest_file')
@click.option('--output-dir', default='data/qiime2', help='Directorio de salida QIIME2')
def import_qiime2(manifest_file, output_dir):
    """Importar datos a QIIME2 desde archivo de manifiesto usando API QIIME2"""
    if not check_qiime2_installation():
        return
    import_to_qiime2(manifest_file, output_dir)


@cli.command()
@click.argument('demux_file')
@click.option('--output-dir', default='results/quality_control', help='Directorio de salida')
@click.option('--min-quality', default=20, help='Calidad mínima para filtrado')
def quality_control(demux_file, output_dir, min_quality):
    """Control de calidad completo: reportes, gráficos y filtrado"""
    qc = QualityControl(demux_file)

    # Ejecutar control de calidad completo
    results = qc.run_quality_control(output_dir, min_quality)

    if results and results['filtered_seqs']:
        print(f"Las secuencias filtradas están listas en: {results['filtered_seqs']}")

@cli.command()
@click.argument('demux_file')
@click.option('--output-dir', default='results/deblur', help='Directorio de salida')
@click.option('--trim-length', default=250, help='Longitud de trim para Deblur')
@click.option('--min-reads', default=10, help='Mínimo de lecturas por muestra para Deblur')
@click.option('--min-size', default=2, help='Mínimo de tamaño para Deblur')
def run_deblur(demux_file, output_dir, trim_length, min_reads, min_size):
    """Ejecutar Deblur para denoising y ASV picking"""
    denoiser = Denoiser(demux_file)

    # Ejecutar Deblur con los parámetros proporcionados
    result = denoiser.run_deblur(
        output_dir=output_dir,
        trim_length=trim_length,
        min_reads=min_reads,
        min_size=min_size
    )

    if result:
        print("🎉 Proceso de denoising con Deblur completado exitosamente!")


@cli.command()
@click.argument('demux_file')
@click.option('--output-dir', default='results/dada2', help='Directorio de salida')
@click.option('--trunc-len-fwd', default=240, help='Punto de corte para forward reads')
@click.option('--trunc-len-rev', default=200, help='Punto de corte para reverse reads')
@click.option('--trim-left-fwd', default=10, help='Bases a remover al inicio (forward)')
@click.option('--trim-left-rev', default=10, help='Bases a remover al inicio (reverse)')
@click.option('--max-ee', default=2.0, help='Expected errors máximo')
def run_dada2(demux_file, output_dir, trunc_len_fwd, trunc_len_rev, trim_left_fwd, trim_left_rev, max_ee):
    """Ejecutar DADA2 para denoising y ASV picking"""
    denoiser = Denoiser(demux_file)

    # Ejecutar DADA2 con los parámetros proporcionados
    result = denoiser.run_dada2(
        output_dir=output_dir,
        trunc_len_fwd=trunc_len_fwd,
        trunc_len_rev=trunc_len_rev,
        trim_left_fwd=trim_left_fwd,
        trim_left_rev=trim_left_rev,
        max_ee_fwd=max_ee,
        max_ee_rev=max_ee
    )

    if result:
        print("🎉 Proceso de denoising con DADA2 completado exitosamente!")


@cli.command()
@click.argument('demux_file')
@click.option('--output-dir', default='results', help='Directorio de salida')
@click.option('--method', default='deblur', help='Método de denoising: deblur o dada2')
def full_analysis(demux_file, output_dir, method):
    """Proceso completo: control de calidad + denoising"""
    # 1. Control de calidad completo
    qc = QualityControl(demux_file)
    print("🎯 Paso 1: Ejecutando control de calidad completo...")
    results = qc.run_quality_control(f"{output_dir}/quality_control")

    # 2. Denoising con las secuencias filtradas
    if results and results['filtered_seqs']:
        print(f"🧬 Paso 2: Ejecutando {method.upper()} con secuencias filtradas...")
        denoiser = Denoiser(results['filtered_seqs'])

        if method == 'deblur':
            denoiser.run_deblur(f"{output_dir}/deblur")
        elif method == 'dada2':
            denoiser.run_dada2(f"{output_dir}/dada2")
        else:
            print(f"❌ Método {method} no reconocido. Usa 'deblur' o 'dada2'")
    else:
        print("❌ No se pudieron obtener secuencias filtradas para el denoising")


if __name__ == '__main__':
    cli()