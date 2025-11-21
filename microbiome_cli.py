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
from modules.taxa import import_database_to_qiime2, taxa_assigner
from modules.phylogeny import make_phylogeny

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
        click.echo("  download               Descargar secuencias SRA desde CSV")
        click.echo("  create-manifest        Crear archivo de manifiesto para QIIME2")
        click.echo("  import-sample-seqs     Importar secuencias a QIIME2")
        click.echo("  quality-control        Control de calidad completo")
        click.echo("  run-deblur             Denoising con Deblur")
        click.echo("  import-reference-database  Importar base de datos de referencia a Qiime2")
        click.echo("  assign-taxonomy        Asignar OTUs/ASVs a taxones")
        click.echo("  build-phylogeny         Generar árbol filogenético")
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
def import_sample_seqs(manifest_file, output_dir):
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

@cli.command()
@click.argument('filename_seq', type=click.Path(exists=True))
@click.argument('filename_taxa', type=click.Path(exists=True))
@click.option('--output-dir', default='ref_database',
              help='Directorio de salida (por defecto: ref_database)')
def import_reference_database(filename_seq, filename_taxa, output_dir):
    """Importar a Qiime2 una base de datos de referencia

    FILENAME_SEQ: Ruta al archivo de secuencias de la base de datos de referencia
    FILENAME_TAXA: Ruta al archivo de taxonomías de la base de datos de referencia

    Ejemplos:
      microbiome_cli.py import-reference-database ref_seqs.fna ref_taxa.txt
      microbiome_cli.py import-reference-database ref_seqs.fna ref_taxa.txt --output-dir my_ref_db
    """
    click.echo(f"📚 Importando base de datos de referencia...")
    click.echo(f"🧬 Secuencias: {filename_seq}")
    click.echo(f"📊 Taxonomías: {filename_taxa}")
    click.echo(f"📁 Directorio de salida: {output_dir}")

    try:
        seq_path, taxa_path = import_database_to_qiime2(filename_seq, filename_taxa, output_dir)
        click.echo(f"✅ Base de datos importada exitosamente:")
        click.echo(f"   - Secuencias: {seq_path}")
        click.echo(f"   - Taxonomías: {taxa_path}")
    except Exception as e:
        click.echo(f"❌ Error importando base de datos: {str(e)}")

@cli.command()
@click.argument('table', type=click.Path(exists=True))
@click.argument('rep_seqs', type=click.Path(exists=True))
@click.argument('seqs_ref', type=click.Path(exists=True))
@click.argument('taxa_ref', type=click.Path(exists=True))
@click.argument('metadata_filename', type=click.Path(exists=True))
@click.option('--cpus', default=1, help='Número de CPUs a usar (por defecto: 1)')
@click.option('--output-dir', default='results/taxonomy',
              help='Directorio de salida (por defecto: results/taxonomy)')
def assign_taxonomy(table, rep_seqs, seqs_ref, taxa_ref, metadata_filename, cpus, output_dir):
    """Asignación taxonómica y generación de archivos CSV por nivel taxonómico

    TABLE: Ruta al artefacto QIIME2 de la tabla de características (.qza)
    REP_SEQS: Ruta al artefacto QIIME2 de secuencias representativas (.qza)
    SEQS_REF: Ruta al artefacto QIIME2 de secuencias de referencia (.qza)
    TAXA_REF: Ruta al artefacto QIIME2 de taxonomía de referencia (.qza)
    METADATA_FILENAME: Ruta al archivo de metadatos (TSV) para QIIME2

    Ejemplos:
      microbiome_cli.py assign-taxonomy table.qza rep-seqs.qza ref-seqs.qza ref-taxa.qza metadata.tsv
      microbiome_cli.py assign-taxonomy table.qza rep-seqs.qza ref-seqs.qza ref-taxa.qza metadata.tsv --cpus 4 --output-dir my_taxa
    """
    click.echo(f"🔍 Asignando taxonomía...")
    click.echo(f"📊 Tabla de características: {table}")
    click.echo(f"🧬 Secuencias representativas: {rep_seqs}")
    click.echo(f"📚 Secuencias de referencia: {seqs_ref}")
    click.echo(f"📚 Taxonomía de referencia: {taxa_ref}")
    click.echo(f"📋 Metadatos: {metadata_filename}")
    click.echo(f"⚡ CPUs: {cpus}")
    click.echo(f"📁 Directorio de salida: {output_dir}")

    try:
        result = taxa_assigner(table, rep_seqs, seqs_ref, taxa_ref, metadata_filename, cpus, output_dir)
        click.echo(f"✅ {result}")
        click.echo(f"📈 Archivos CSV generados:")
        click.echo(f"   - phylum.csv, class.csv, order.csv")
        click.echo(f"   - family.csv, genus.csv, species.csv")
        click.echo(f"   - taxa_barplot.qzv (visualización QIIME2)")
    except Exception as e:
        click.echo(f"❌ Error en asignación taxonómica: {str(e)}")

@cli.command()
@click.argument('rep_seqs', type=click.Path(exists=True))
@click.option('--output-dir', default='results/phylogeny',
              help='Directorio de salida (por defecto: results/phylogeny)')
def build_phylogeny(rep_seqs, output_dir):
    """Generar árbol filogenético a partir de secuencias representativas

    REP_SEQS: Ruta al artefacto QIIME2 de secuencias representativas (.qza)

    Genera dos árboles filogenéticos:
      - unrooted_tree.qza: Árbol filogenético sin raíz
      - rooted_tree.qza: Árbol filogenético con raíz

    El proceso incluye:
      - Alineamiento múltiple con MAFFT
      - Construcción de árbol con FastTree
      - Enraizamiento del árbol

    Ejemplos:
      microbiome_cli.py make-phylogeny rep-seqs.qza
      microbiome_cli.py make-phylogeny rep-seqs.qza --output-dir my_phylogeny
    """
    click.echo(f"🌳 Generando árbol filogenético...")
    click.echo(f"🧬 Secuencias representativas: {rep_seqs}")
    click.echo(f"📁 Directorio de salida: {output_dir}")

    try:
        unrooted_path, rooted_path = make_phylogeny(rep_seqs, output_dir)
        click.echo(f"✅ Árbol filogenético generado exitosamente:")
        click.echo(f"   - Árbol sin raíz: {unrooted_path}")
        click.echo(f"   - Árbol con raíz: {rooted_path}")
        click.echo(f"🌿 Los árboles están listos para análisis de diversidad filogenética")
    except Exception as e:
        click.echo(f"❌ Error generando árbol filogenético: {str(e)}")



if __name__ == '__main__':
    cli()