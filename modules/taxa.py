# modules/taxa.py
import pandas as pd
import tempfile
import pathlib
import os
from qiime2 import Artifact, Metadata
from qiime2.plugins.feature_classifier.pipelines import classify_consensus_vsearch
from qiime2.plugins.taxa.visualizers import barplot


def import_database_to_qiime2(filename_seqs, filename_taxa, output_dir):
    """Importar base de datos de referencia a Qiime2"""
    os.makedirs(output_dir, exist_ok=True)

    seq = Artifact.import_data('FeatureData[Sequence]', filename_seqs)
    taxa = Artifact.import_data('FeatureData[Taxonomy]', filename_taxa, 'HeaderlessTSVTaxonomyFormat')

    seq.save(f"{output_dir}/reference_sequences.qza")
    taxa.save(f"{output_dir}/reference_taxonomy.qza")

    return f"{output_dir}/reference_sequences.qza", f"{output_dir}/reference_taxonomy.qza"


def load_reference_DB_artifact(filename_seqs_artifact, filename_taxa_artifact):
    """Cargar artefactos de base de datos de referencia"""
    seq = Artifact.load(filename_seqs_artifact)
    taxa = Artifact.load(filename_taxa_artifact)
    return seq, taxa


def normalized_df(dataframe):
    """Normalizar dataframe a porcentajes"""
    columns = dataframe.columns
    for column in columns:
        dataframe[column] = (dataframe[column] / dataframe[column].sum()) * 100
    return dataframe


def taxa_assigner(table, rep_seqs, seqs_ref, taxa_ref, metadata_filename, cpus, output_folder):
    """Asignar taxonomía y generar archivos CSV por nivel taxonómico"""
    os.makedirs(output_folder, exist_ok=True)

    # Cargar artefactos si se pasan como paths
    if isinstance(table, str):
        table = Artifact.load(table)
    if isinstance(rep_seqs, str):
        rep_seqs = Artifact.load(rep_seqs)
    if isinstance(seqs_ref, str):
        seqs_ref = Artifact.load(seqs_ref)
    if isinstance(taxa_ref, str):
        taxa_ref = Artifact.load(taxa_ref)

    # Clasificación taxonómica
    taxonomy = classify_consensus_vsearch(
        query=rep_seqs,
        reference_reads=seqs_ref,
        reference_taxonomy=taxa_ref,
        threads=cpus
    )

    # Crear barplot de taxonomía
    taxa_barplot = barplot(
        table=table,
        taxonomy=taxonomy.classification,
        metadata=Metadata.load(metadata_filename)
    )
    taxa_barplot = taxa_barplot.visualization
    taxa_barplot.save(f"{output_folder}/taxa_barplot.qzv")

    # Exportar datos y generar archivos CSV por nivel taxonómico
    csvs_barplot = []
    with tempfile.TemporaryDirectory() as tmpdir:
        taxa_barplot.export_data(tmpdir)
        data_dir_fp = pathlib.Path(tmpdir)
        csv_fps = sorted(data_dir_fp.glob('level-*csv'))

        for csv_fp in csv_fps:
            df_barplot = pd.read_csv(csv_fp, index_col='index')
            csvs_barplot.append(df_barplot)

    # Generar archivos CSV para cada nivel taxonómico
    levels = {
        1: "phylum",
        2: "class",
        3: "order",
        4: "family",
        5: "genus",
        -1: "species"  # Último nivel para especies
    }

    for level_idx, level_name in levels.items():
        df_level = csvs_barplot[level_idx].T
        df_level.drop(df_level.index[-1], inplace=True)  # Eliminar fila "Unknown"
        df_level = normalized_df(df_level)
        df_level.to_csv(f"{output_folder}/{level_name}.csv")

    return f"Archivos taxonómicos generados en: {output_folder}"