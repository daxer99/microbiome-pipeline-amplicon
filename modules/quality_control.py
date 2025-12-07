"""
Módulo simplificado para control de calidad - Solo gráficos de calidad
"""
import warnings

warnings.filterwarnings("ignore", category=UserWarning)  # Silenciar warnings específicos

import matplotlib

matplotlib.use('Agg')  # Usar backend no interactivo
import matplotlib.pyplot as plt
from pathlib import Path
from qiime2 import Artifact
from qiime2.plugins.demux.visualizers import summarize
import dokdo
import tempfile
import os


class QualityControl:
    """Clase simplificada para control de calidad con dokdo"""

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

    def run_quality_control(self, output_dir="results/quality_control"):
        """
        Ejecuta el control de calidad simplificado:
        1. Genera visualización de calidad QIIME2
        2. Crea gráficos de perfil de calidad con dokdo

        Args:
            output_dir: Directorio de salida
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print("🎯 Iniciando control de calidad...")

        # Paso 1: Visualización de calidad QIIME2
        print("📊 Generando visualización de calidad QIIME2...")
        viz_path = self.create_quality_visualization(output_dir)

        # Paso 2: Gráficos de calidad con dokdo
        print("📈 Creando gráficos de perfil de calidad...")
        plot_path = self.plot_quality_profile(output_path / "quality_profile.png")

        print("✅ Control de calidad completado exitosamente!")

        return {
            'quality_viz': viz_path,
            'quality_plot': plot_path
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

    def plot_quality_profile(self, output_file="quality_profile.png"):
        """Genera gráficos de perfil de calidad usando dokdo"""
        if self.quality_visualization is None:
            self.create_quality_visualization()

        print(f"📊 Generando gráfico de calidad...")

        viz = self.quality_visualization.visualization

        # Verificar si es paired-end
        is_paired = False
        try:
            # Intentar crear un gráfico con ambos strands
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            # Configurar warnings temporalmente
            import warnings
            from matplotlib.cbook import mplDeprecation
            warnings.filterwarnings("ignore", category=UserWarning)
            warnings.filterwarnings("ignore", category=mplDeprecation)

            # Intentar forward
            try:
                dokdo.read_quality_plot(viz, strand='forward', ax=axes[0])
                axes[0].set_title('Forward Read', fontsize=12, fontweight='bold')
                axes[0].set_xlabel('Position in Read', fontsize=10)
                axes[0].set_ylabel('Quality Score', fontsize=10)
            except Exception:
                axes[0].set_visible(False)

            # Intentar reverse
            try:
                dokdo.read_quality_plot(viz, strand='reverse', ax=axes[1])
                axes[1].set_title('Reverse Read', fontsize=12, fontweight='bold')
                axes[1].set_xlabel('Position in Read', fontsize=10)
                axes[1].set_ylabel('')
                is_paired = True
            except Exception:
                axes[1].set_visible(False)
                is_paired = False

            # Si solo hay un gráfico visible, ajustar
            if not axes[1].get_visible() and axes[0].get_visible():
                # Single-end: usar figura más pequeña
                plt.close(fig)
                fig, ax = plt.subplots(1, 1, figsize=(10, 6))
                dokdo.read_quality_plot(viz, strand='forward', ax=ax)
                ax.set_title('Quality Profile', fontsize=14, fontweight='bold')
                ax.set_xlabel('Position in Read', fontsize=11)
                ax.set_ylabel('Quality Score', fontsize=11)

        except Exception as e:
            print(f"⚠️  Advertencia: {e}")
            # Fallback: gráfico simple
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            dokdo.read_quality_plot(viz, strand='forward', ax=ax)
            ax.set_title('Quality Profile', fontsize=14, fontweight='bold')
            ax.set_xlabel('Position in Read', fontsize=11)
            ax.set_ylabel('Quality Score', fontsize=11)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"✅ Gráfico de calidad guardado: {output_file}")

        return output_file

    def get_visualization(self):
        """Retorna la visualización de calidad"""
        return self.quality_visualization