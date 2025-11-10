"""
Módulo para control de calidad - Gráficos y filtrado
"""
import matplotlib.pyplot as plt
from pathlib import Path
from qiime2 import Artifact
from qiime2.plugins.demux.visualizers import summarize

try:
    import dokdo

    DOKDO_AVAILABLE = True
except ImportError:
    DOKDO_AVAILABLE = False
    print("⚠️  dokdo no está instalado. Los gráficos de calidad no estarán disponibles.")
    print("   Instala con: pip install dokdo")


class QualityControl:
    """Clase para control de calidad completo"""

    def __init__(self, demux_artifact):
        """
        Args:
            demux_artifact: Path al archivo .qza o objeto Artifact de QIIME2
        """
        if isinstance(demux_artifact, (str, Path)):
            self.demux_seqs = Artifact.load(str(demux_artifact))
        else:
            self.demux_seqs = demux_artifact

        self.quality_visualization = None
        self.filtered_seqs = None

    def run_quality_control(self, output_dir="results/quality_control", min_quality=20):
        """
        Ejecuta el control de calidad completo:
        1. Genera visualización de calidad
        2. Crea gráficos de perfil de calidad
        3. Aplica filtrado de calidad

        Args:
            output_dir: Directorio de salida
            min_quality: Calidad mínima para filtrado
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print("🎯 Iniciando control de calidad completo...")

        # Paso 1: Visualización de calidad QIIME2
        print("📊 Paso 1/3: Generando visualización de calidad QIIME2...")
        viz_path = self.create_quality_visualization(output_dir)

        # Paso 2: Gráficos de calidad
        print("📈 Paso 2/3: Creando gráficos de perfil de calidad...")
        plot_path = self.plot_quality_profile(output_path / "quality_profile.png")

        # Paso 3: Filtrado de calidad
        print("🔧 Paso 3/3: Aplicando filtrado de calidad...")
        filtered_result = self.run_quality_filter(output_dir, min_quality)

        if filtered_result:
            self.filtered_seqs = Artifact.load(str(filtered_result['filtered_seqs']))

        print("✅ Control de calidad completado exitosamente!")

        return {
            'quality_viz': viz_path,
            'quality_plot': plot_path,
            'filtered_seqs': filtered_result['filtered_seqs'] if filtered_result else None,
            'filter_stats': filtered_result['filter_stats'] if filtered_result else None
        }

    def create_quality_visualization(self, output_dir="results/quality"):
        """Crea la visualización de calidad de QIIME2"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print("🔍 Generando visualización de calidad...")
        self.quality_visualization = summarize(self.demux_seqs)

        # Guardar visualización
        viz_path = output_path / "quality_viz.qzv"
        self.quality_visualization.visualization.save(str(viz_path))
        print(f"✅ Visualización de calidad guardada: {viz_path}")

        return viz_path

    def plot_quality_profile(self, output_file="quality_profile.png", figsize=(15, 6)):
        """Genera gráficos de perfil de calidad usando dokdo"""
        if not DOKDO_AVAILABLE:
            print("❌ dokdo no está disponible para generar gráficos")
            return None

        if self.quality_visualization is None:
            self.create_quality_visualization()

        print(f"📊 Generando gráfico de calidad...")

        fig = plt.figure(figsize=figsize)

        # Asumimos paired-end por ahora, pero podríamos detectarlo
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)

        dokdo.read_quality_plot(self.quality_visualization.visualization,
                                strand='forward', ax=ax1)
        dokdo.read_quality_plot(self.quality_visualization.visualization,
                                strand='reverse', ax=ax2)

        ax1.set_title('Read Forward - Calidad', fontsize=12, fontweight='bold')
        ax2.set_title('Read Reverse - Calidad', fontsize=12, fontweight='bold')

        # Ajustes estéticos
        for ax in [ax1, ax2]:
            ax.tick_params(axis='both', which='major', labelsize=10)
            ax.set_xlabel('Posición en el read', fontsize=11)
            ax.set_ylabel('Score de Calidad', fontsize=11)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✅ Gráfico de calidad guardado: {output_file}")

        return output_file

    def run_quality_filter(self, output_dir="results/quality_filtered", min_quality=20):
        """
        Filtrado de calidad básico

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

    def get_filtered_seqs(self):
        """Retorna las secuencias filtradas"""
        return self.filtered_seqs