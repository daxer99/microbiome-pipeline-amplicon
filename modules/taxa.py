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

    print("🔍 DEBUG: Iniciando carga de artefactos...")
    # Cargar artefactos si se pasan como paths
    if isinstance(table, str):
        print(f"   Cargando table desde: {table}")
        table = Artifact.load(table)
    if isinstance(rep_seqs, str):
        print(f"   Cargando rep_seqs desde: {rep_seqs}")
        rep_seqs = Artifact.load(rep_seqs)
    if isinstance(seqs_ref, str):
        print(f"   Cargando seqs_ref desde: {seqs_ref}")
        seqs_ref = Artifact.load(seqs_ref)
    if isinstance(taxa_ref, str):
        print(f"   Cargando taxa_ref desde: {taxa_ref}")
        taxa_ref = Artifact.load(taxa_ref)
    print("✅ DEBUG: Artefactos cargados correctamente")

    # Clasificación taxonómica
    print("🔍 DEBUG: Iniciando clasificación taxonómica con VSEARCH...")
    try:
        taxonomy = classify_consensus_vsearch(
            query=rep_seqs,
            reference_reads=seqs_ref,
            reference_taxonomy=taxa_ref,
            threads=cpus
        )
        print("✅ DEBUG: Clasificación taxonómica completada")
    except Exception as e:
        print(f"❌ DEBUG: Error en classify_consensus_vsearch: {e}")
        import traceback
        print(traceback.format_exc())
        raise

    # Guardar la clasificación taxonómica
    print(f"🔍 DEBUG: Guardando clasificación en: {output_folder}/taxonomy.qza")
    try:
        taxonomy.classification.save(f"{output_folder}/taxonomy.qza")
        print("✅ DEBUG: Clasificación guardada correctamente")
    except Exception as e:
        print(f"❌ DEBUG: Error al guardar clasificación: {e}")
        import traceback
        print(traceback.format_exc())
        raise

    # Crear barplot de taxonomía
    print("🔍 DEBUG: Creando barplot de taxonomía...")
    try:
        taxa_barplot = barplot(
            table=table,
            taxonomy=taxonomy.classification,
            metadata=Metadata.load(metadata_filename)
        )
        taxa_barplot = taxa_barplot.visualization
        print("✅ DEBUG: Barplot creado correctamente")
    except Exception as e:
        print(f"❌ DEBUG: Error al crear barplot: {e}")
        import traceback
        print(traceback.format_exc())
        raise

    print(f"🔍 DEBUG: Guardando barplot en: {output_folder}/taxa_barplot.qzv")
    try:
        taxa_barplot.save(f"{output_folder}/taxa_barplot.qzv")
        print("✅ DEBUG: Barplot guardado correctamente")
    except Exception as e:
        print(f"❌ DEBUG: Error al guardar barplot: {e}")
        import traceback
        print(traceback.format_exc())
        raise

    # Exportar datos y generar archivos CSV por nivel taxonómico
    print("🔍 DEBUG: Exportando datos de barplot...")
    csvs_barplot = []
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            taxa_barplot.export_data(tmpdir)
            print(f"✅ DEBUG: Datos exportados a: {tmpdir}")

            data_dir_fp = pathlib.Path(tmpdir)
            csv_fps = sorted(data_dir_fp.glob('level-*.csv'))
            print(f"🔍 DEBUG: Archivos CSV encontrados: {[str(fp) for fp in csv_fps]}")

            for csv_fp in csv_fps:
                df_barplot = pd.read_csv(csv_fp, index_col='index')
                csvs_barplot.append(df_barplot)
                print(f"✅ DEBUG: CSV cargado: {csv_fp.name} - Shape: {df_barplot.shape}")
        except Exception as e:
            print(f"❌ DEBUG: Error al exportar datos: {e}")
            import traceback
            print(traceback.format_exc())
            raise

    print(f"🔍 DEBUG: Total de niveles taxonómicos encontrados: {len(csvs_barplot)}")

    # Definir niveles taxonómicos según QIIME2
    levels = {
        1: "phylum",
        2: "class",
        3: "order",
        4: "family",
        5: "genus",
        6: "species"
    }

    # Generar archivos CSV para cada nivel taxonómico
    print("🔍 DEBUG: Generando archivos CSV por nivel taxonómico...")
    for level_idx, level_name in levels.items():
        print(f"   Procesando nivel {level_idx} ({level_name})...")
        # Verificar que el índice existe en csvs_barplot
        if level_idx < len(csvs_barplot):
            try:
                df_level = csvs_barplot[level_idx].T
                print(f"      DataFrame transpuesto - Shape: {df_level.shape}")

                # Eliminar fila "Unknown" si existe
                if len(df_level.index) > 0 and df_level.index[-1] == "Unknown":
                    df_level = df_level.drop(df_level.index[-1])
                    print(f"      Fila 'Unknown' eliminada")

                # Normalizar a porcentajes
                df_level = normalized_df(df_level)
                print(f"      DataFrame normalizado")

                # Guardar CSV
                output_file = f"{output_folder}/{level_name}.csv"
                df_level.to_csv(output_file)
                print(f"✅ Archivo generado: {output_file}")
            except Exception as e:
                print(f"❌ DEBUG: Error procesando nivel {level_name}: {e}")
                import traceback
                print(traceback.format_exc())
                raise
        else:
            print(f"⚠️  Advertencia: No se encontró nivel {level_name} (índice {level_idx}) en los resultados")

    return f"Archivos taxonómicos generados en: {output_folder}"