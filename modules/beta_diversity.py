# modules/beta_diversity.py
import tempfile
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import os
from qiime2 import Artifact
from qiime2.plugins.diversity.pipelines import beta, beta_phylogenetic
from qiime2.plugins.diversity.methods import pcoa
import dokdo


def calculate_beta_diversity(table, metrics, output_dir):
    """Calcula matrices de distancia beta no filogenéticas.

    Args:
        table: Ruta al artefacto QIIME2 de la tabla de características o el artefacto mismo.
        metrics: Lista de métricas de distancia beta.
        output_dir: Directorio donde guardar los resultados.

    Returns:
        list: Rutas a las matrices de distancia guardadas
    """
    os.makedirs(output_dir, exist_ok=True)

    if isinstance(table, str):
        table = Artifact.load(table)

    distance_matrices = []
    output_files = []

    for metric in metrics:
        beta_calculator = beta(table=table, metric=metric, n_jobs='auto')
        distance_matrix = beta_calculator.distance_matrix

        # Guardar matriz de distancia
        matrix_path = f"{output_dir}/{metric}_distance_matrix.qza"
        distance_matrix.save(matrix_path)
        distance_matrices.append(distance_matrix)

        # Exportar a CSV
        with tempfile.TemporaryDirectory() as tmpdir:
            distance_matrix.export_data(tmpdir)
            beta_dir_fp = pathlib.Path(tmpdir)
            csv_file = list(beta_dir_fp.glob('*.tsv'))[0]
            df = pd.read_table(csv_file, index_col=0)
            csv_path = f"{output_dir}/{metric}_distance_matrix.csv"
            df.to_csv(csv_path)
            output_files.append(csv_path)

    return output_files, distance_matrices


def calculate_phylogenetic_beta_diversity(table, metrics, rooted_tree, output_dir):
    """Calcula matrices de distancia beta filogenéticas.

    Args:
        table: Ruta al artefacto QIIME2 de la tabla de características o el artefacto mismo.
        metrics: Lista de métricas de distancia beta filogenéticas.
        rooted_tree: Ruta al artefacto QIIME2 del árbol filogenético enraizado.
        output_dir: Directorio donde guardar los resultados.

    Returns:
        list: Rutas a las matrices de distancia guardadas
    """
    os.makedirs(output_dir, exist_ok=True)

    if isinstance(table, str):
        table = Artifact.load(table)
    if isinstance(rooted_tree, str):
        rooted_tree = Artifact.load(rooted_tree)

    distance_matrices = []
    output_files = []

    for metric in metrics:
        beta_calculator = beta_phylogenetic(
            table=table,
            metric=metric,
            threads='auto',
            phylogeny=rooted_tree
        )
        distance_matrix = beta_calculator.distance_matrix

        # Guardar matriz de distancia
        matrix_path = f"{output_dir}/{metric}_distance_matrix.qza"
        distance_matrix.save(matrix_path)
        distance_matrices.append(distance_matrix)

        # Exportar a CSV
        with tempfile.TemporaryDirectory() as tmpdir:
            distance_matrix.export_data(tmpdir)
            beta_dir_fp = pathlib.Path(tmpdir)
            csv_file = list(beta_dir_fp.glob('*.tsv'))[0]
            df = pd.read_table(csv_file, index_col=0)
            csv_path = f"{output_dir}/{metric}_distance_matrix.csv"
            df.to_csv(csv_path)
            output_files.append(csv_path)

    return output_files, distance_matrices


def plot_pcoa(distance_matrix, metadata, hue, output_file, metric_name):
    """Genera gráfico PCoA a partir de una matriz de distancia.

    Args:
        distance_matrix: Artefacto QIIME2 de matriz de distancia.
        metadata: Ruta al archivo de metadatos.
        hue: Columna en los metadatos para colorear los puntos.
        output_file: Ruta del archivo de salida para el gráfico.
        metric_name: Nombre de la métrica de distancia.
    """
    if isinstance(distance_matrix, str):
        distance_matrix = Artifact.load(distance_matrix)

    # Realizar PCoA usando el método de QIIME2
    pcoa_result = pcoa(distance_matrix)

    # Crear gráfico usando dokdo
    fig, ax = plt.subplots(figsize=(8, 8))
    dokdo.beta_2d_plot(pcoa_result.pcoa, metadata=metadata, hue=hue, ax=ax)
    ax.set_title(f'PCoA - {metric_name}')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    return output_file