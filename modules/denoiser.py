"""
Módulo para denoising con DADA2 y Deblur
"""
from pathlib import Path
from qiime2 import Artifact


class Denoiser:
    """Clase para denoising con DADA2 y Deblur"""

    def __init__(self, demux_artifact):
        """
        Args:
            demux_artifact: Path al archivo .qza o objeto Artifact de QIIME2
        """
        if isinstance(demux_artifact, (str, Path)):
            self.demux_seqs = Artifact.load(str(demux_artifact))
        else:
            self.demux_seqs = demux_artifact

    def run_dada2(self, output_dir="results/dada2", **dada2_params):
        """
        Ejecuta DADA2 para denoising y obtención de ASVs

        Args:
            output_dir: Directorio de salida
            dada2_params: Parámetros para DADA2
        """
        from qiime2.plugins.dada2.methods import denoise_paired, denoise_single

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Parámetros por defecto
        default_params = {
            'trunc_len_fwd': 240,
            'trunc_len_rev': 200,
            'trim_left_fwd': 10,
            'trim_left_rev': 10,
            'max_ee_fwd': 2.0,
            'max_ee_rev': 2.0,
            'trunc_q': 2,
            'min_overlap': 12,
            'n_threads': 4
        }

        # Actualizar con parámetros proporcionados
        default_params.update(dada2_params)

        print("🧬 Ejecutando DADA2 para denoising...")
        print("📋 Parámetros utilizados:")
        for key, value in default_params.items():
            print(f"   • {key}: {value}")

        try:
            # Por simplicidad, asumimos paired-end
            dada2_result = denoise_paired(
                demultiplexed_seqs=self.demux_seqs,
                trunc_len_fwd=default_params['trunc_len_fwd'],
                trunc_len_rev=default_params['trunc_len_rev'],
                trim_left_fwd=default_params['trim_left_fwd'],
                trim_left_rev=default_params['trim_left_rev'],
                max_ee_fwd=default_params['max_ee_fwd'],
                max_ee_rev=default_params['max_ee_rev'],
                trunc_q=default_params['trunc_q'],
                min_overlap=default_params['min_overlap'],
                n_threads=default_params['n_threads']
            )

            # Guardar resultados
            table = dada2_result.table
            rep_seqs = dada2_result.representative_sequences
            stats = dada2_result.denoising_stats

            table_path = output_path / "table.qza"
            rep_seqs_path = output_path / "rep-seqs.qza"
            stats_path = output_path / "stats.qza"

            table.save(str(table_path))
            rep_seqs.save(str(rep_seqs_path))
            stats.save(str(stats_path))

            print("✅ Denoising con DADA2 completado:")
            print(f"   • Tabla de features: {table_path}")
            print(f"   • Secuencias representativas: {rep_seqs_path}")
            print(f"   • Estadísticas: {stats_path}")

            return {
                'table': table_path,
                'rep_seqs': rep_seqs_path,
                'stats': stats_path
            }

        except Exception as e:
            print(f"❌ Error en DADA2: {e}")
            return None

    def run_deblur(self, output_dir="results/deblur", **deblur_params):
        """
        Ejecuta Deblur para denoising y obtención de ASVs

        Args:
            output_dir: Directorio de salida
            deblur_params: Parámetros para Deblur
        """
        from qiime2.plugins.deblur.methods import denoise_16S

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Parámetros por defecto para Deblur
        default_params = {
            'trim_length': 250,
            'min_reads': 10,
            'min_size': 2,
            'jobs_to_start': 4,
            'hashed_feature_ids': True
        }

        # Actualizar con parámetros proporcionados
        default_params.update(deblur_params)

        print("🧹 Ejecutando Deblur para denoising...")
        print("📋 Parámetros utilizados:")
        for key, value in default_params.items():
            print(f"   • {key}: {value}")

        try:
            # Ejecutar Deblur
            deblur_result = denoise_16S(
                demultiplexed_seqs=self.demux_seqs,
                trim_length=default_params['trim_length'],
                min_reads=default_params['min_reads'],
                min_size=default_params['min_size'],
                jobs_to_start=default_params['jobs_to_start'],
                hashed_feature_ids=default_params['hashed_feature_ids']
            )

            # Guardar resultados
            table = deblur_result.table
            rep_seqs = deblur_result.representative_sequences
            stats = deblur_result.deblur_stats

            table_path = output_path / "table.qza"
            rep_seqs_path = output_path / "rep-seqs.qza"
            stats_path = output_path / "stats.qza"

            table.save(str(table_path))
            rep_seqs.save(str(rep_seqs_path))
            stats.save(str(stats_path))

            print("✅ Denoising con Deblur completado:")
            print(f"   • Tabla de features: {table_path}")
            print(f"   • Secuencias representativas: {rep_seqs_path}")
            print(f"   • Estadísticas: {stats_path}")

            return {
                'table': table_path,
                'rep_seqs': rep_seqs_path,
                'stats': stats_path
            }

        except Exception as e:
            print(f"❌ Error en Deblur: {e}")
            return None

    def run_quality_filter(self, output_dir="results/quality_filtered", min_quality=20):
        """
        Filtrado de calidad básico antes de Deblur o DADA2

        Args:
            output_dir: Directorio de salida
            min_quality: Calidad mínima promedio
        """
        from qiime2.plugins.quality_filter.methods import q_score

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print("🔧 Aplicando filtrado de calidad...")

        try:
            # Filtrado por calidad
            quality_result = q_score(
                demux=self.demux_seqs,
                min_quality=min_quality,
                quality_window=3,
                min_length_fraction=0.75,
                max_ambiguous=0
            )

            # Guardar resultados filtrados
            filtered_seqs = quality_result.filtered_sequences
            filtered_stats = quality_result.filter_stats

            filtered_path = output_path / "filtered_seqs.qza"
            stats_path = output_path / "filter_stats.qza"

            filtered_seqs.save(str(filtered_path))
            filtered_stats.save(str(stats_path))

            print("✅ Filtrado de calidad completado:")
            print(f"   • Secuencias filtradas: {filtered_path}")
            print(f"   • Estadísticas de filtrado: {stats_path}")

            return {
                'filtered_seqs': filtered_path,
                'filter_stats': stats_path
            }

        except Exception as e:
            print(f"❌ Error en filtrado de calidad: {e}")
            return None