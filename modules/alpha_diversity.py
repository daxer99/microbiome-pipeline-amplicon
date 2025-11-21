# modules/alpha_diversity.py
import tempfile
import pathlib
import pandas as pd
import os
from qiime2 import Artifact
from qiime2.plugins.diversity.pipelines import alpha, alpha_phylogenetic


def calculate_alpha_diversity(table, metrics, output_folder, rooted_tree=None):
    """Calcula la diversidad alfa para una lista de métricas.

    Args:
        table: Ruta al artefacto QIIME2 de la tabla de características o el artefacto mismo.
        metrics: Lista de métricas de diversidad alfa a calcular.
        output_folder: Directorio donde se guardarán los archivos CSV resultantes.
        rooted_tree: Ruta al artefacto QIIME2 del árbol filogenético enraizado (necesario para faith_pd).

    Returns:
        list: Rutas a los archivos CSV generados
    """
    os.makedirs(output_folder, exist_ok=True)

    # Cargar artefactos si se pasan como strings
    if isinstance(table, str):
        table = Artifact.load(table)
    if rooted_tree and isinstance(rooted_tree, str):
        rooted_tree = Artifact.load(rooted_tree)

    diversity = []
    for metric in metrics:
        if metric == 'faith_pd':
            if not rooted_tree:
                raise ValueError("La métrica 'faith_pd' requiere un árbol filogenético enraizado.")
            alpha_calculator = alpha_phylogenetic(table=table, phylogeny=rooted_tree, metric='faith_pd')
        else:
            alpha_calculator = alpha(table=table, metric=metric)

        with tempfile.TemporaryDirectory() as tmpdir:
            alpha_calculator.alpha_diversity.export_data(tmpdir)
            alpha_dir_fp = pathlib.Path(tmpdir)
            csv_file = sorted(alpha_dir_fp.glob('*.tsv'))
            for csv in csv_file:
                df_lista = pd.read_table(csv, index_col=False)
                diversity.append(df_lista)

    # Guardar cada métrica en un archivo CSV
    output_files = []
    for i, metric in enumerate(metrics):
        data = diversity[i]
        diversity_path = f"{output_folder}/{metric}.csv"
        data.to_csv(diversity_path, index=False)
        output_files.append(diversity_path)

    return output_files