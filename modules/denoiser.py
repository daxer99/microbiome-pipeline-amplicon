"""
Módulo para denoising con Deblur
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
            'left_trim_len':0,
            'trim_length': 250,
            'min_reads': 50,
            'min_size': 2,
            'jobs_to_start': 8,
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
                left_trim_len =default_params['left_trim_len'],
                trim_length=default_params['trim_length'],
                min_reads=default_params['min_reads'],
                min_size=default_params['min_size'],
                jobs_to_start=default_params['jobs_to_start'],
                sample_stats=True,
            )

            # Guardar resultados
            table = deblur_result.table
            rep_seqs = deblur_result.representative_sequences
            stats = deblur_result.stats

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