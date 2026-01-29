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

    # Guardar la clasificación taxonómica
    taxonomy.classification.save(f"{output_folder}/taxonomy.qza")

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
        csv_fps = sorted(data_dir_fp.glob('level-*.csv'))

        for csv_fp in csv_fps:
            df_barplot = pd.read_csv(csv_fp, index_col='index')
            csvs_barplot.append(df_barplot)

    # Mapeo correcto de índices (QIIME2 empieza en level-1, level-2, etc.)
    # Los archivos CSV están ordenados, así que:
    # csvs_barplot[0] = level-1 (Kingdom/Domain)
    # csvs_barplot[1] = level-2 (Phylum)
    # csvs_barplot[2] = level-3 (Class)
    # csvs_barplot[3] = level-4 (Order)
    # csvs_barplot[4] = level-5 (Family)
    # csvs_barplot[5] = level-6 (Genus)
    # csvs_barplot[6] = level-7 (Species)

    # Definir niveles taxonómicos según QIIME2
    levels = {
        1: "phylum",  # csvs_barplot[1]
        2: "class",  # csvs_barplot[2]
        3: "order",  # csvs_barplot[3]
        4: "family",  # csvs_barplot[4]
        5: "genus",  # csvs_barplot[5]
        6: "species"  # csvs_barplot[6]
    }

    # Generar archivos CSV para cada nivel taxonómico
    for level_idx, level_name in levels.items():
        # Verificar que el índice existe en csvs_barplot
        if level_idx < len(csvs_barplot):
            df_level = csvs_barplot[level_idx].T

            # Eliminar fila "Unknown" si existe
            if df_level.index[-1] == "Unknown":
                df_level = df_level.drop(df_level.index[-1])

            # Normalizar a porcentajes
            df_level = normalized_df(df_level)

            # Guardar CSV
            output_file = f"{output_folder}/{level_name}.csv"
            df_level.to_csv(output_file)
            print(f"✅ Archivo generado: {output_file}")
        else:
            print(f"⚠️  Advertencia: No se encontró nivel {level_name} (índice {level_idx}) en los resultados")

    return f"Archivos taxonómicos generados en: {output_folder}"