"""
Módulo simplificado para control de calidad - Solo gráficos de calidad
"""
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

    def _detect_strands(self):
        """Detecta si es single-end o paired-end"""
        if self.quality_visualization is None:
            self.create_quality_visualization()

        viz = self.quality_visualization.visualization

        # Detectar strands disponibles
        with tempfile.TemporaryDirectory() as tmpdir:
            viz.export_data(tmpdir)

            forward_exists = os.path.exists(os.path.join(tmpdir, 'forward-seven-number-summaries.tsv'))
            reverse_exists = os.path.exists(os.path.join(tmpdir, 'reverse-seven-number-summaries.tsv'))

            if forward_exists and reverse_exists:
                return 'paired'
            elif forward_exists:
                return 'single'
            else:
                return None

    def plot_quality_profile(self, output_file="quality_profile.png"):
        """Genera gráficos de perfil de calidad usando dokdo"""
        if self.quality_visualization is None:
            self.create_quality_visualization()

        print(f"📊 Generando gráfico de calidad...")

        # Detectar tipo de datos
        data_type = self._detect_strands()

        if not data_type:
            print("⚠️  No se pudieron detectar datos de calidad")
            return None

        viz = self.quality_visualization.visualization

        # Configurar figura según el tipo de datos
        if data_type == 'paired':
            # Paired-end: crear figura con 2 subplots
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            # Crear gráfico forward
            ax1 = axes[0]
            dokdo.read_quality_plot(viz, strand='forward', ax=ax1)
            ax1.set_title('Forward Read - Quality Profile', fontsize=12, fontweight='bold')
            ax1.set_xlabel('Position in Read', fontsize=10)
            ax1.set_ylabel('Quality Score', fontsize=10)
            ax1.grid(True, alpha=0.3)

            # Crear gráfico reverse
            ax2 = axes[1]
            dokdo.read_quality_plot(viz, strand='reverse', ax=ax2)
            ax2.set_title('Reverse Read - Quality Profile', fontsize=12, fontweight='bold')
            ax2.set_xlabel('Position in Read', fontsize=10)
            ax2.set_ylabel('')
            ax2.grid(True, alpha=0.3)

            # Ajustar layout
            plt.tight_layout()

        else:  # single-end
            # Single-end: crear figura con 1 subplot
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))

            dokdo.read_quality_plot(viz, strand='forward', ax=ax)
            ax.set_title('Quality Profile', fontsize=14, fontweight='bold')
            ax.set_xlabel('Position in Read', fontsize=11)
            ax.set_ylabel('Quality Score', fontsize=11)
            ax.grid(True, alpha=0.3)

            plt.tight_layout()

        # Guardar figura
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"✅ Gráfico de calidad guardado: {output_file}")

        return output_file

    def get_visualization(self):
        """Retorna la visualización de calidad"""
        return self.quality_visualization