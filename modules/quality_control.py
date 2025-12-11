"""
Módulo para control de calidad con dokdo y filtrado opcional por q-score
"""
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from qiime2 import Artifact, Visualization
from qiime2.plugins.demux.visualizers import summarize
from qiime2.plugins.quality_filter.methods import q_score
import os

# Intentar importar dokdo
try:
    import dokdo

    DOKDO_AVAILABLE = True
except ImportError:
    DOKDO_AVAILABLE = False
    print("⚠️  Warning: dokdo not available. Install with: pip install dokdo")


class QualityControl:
    """Clase para control de calidad usando dokdo"""

    def __init__(self, demux_artifact):
        """
        Args:
            demux_artifact: Path al archivo .qza o objeto Artifact de QIIME2
        """
        if isinstance(demux_artifact, (str, Path)):
            self.demux_seqs = Artifact.load(str(demux_artifact))
            self.demux_path = str(demux_artifact)
        else:
            self.demux_seqs = demux_artifact
            self.demux_path = None

        self.quality_visualization = None
        self.quality_viz_path = None

    def run_quality_control(self, output_dir="results/quality_control",
                            filter_sequences=False, min_quality=4):
        """
        Ejecuta el control de calidad:
        1. Genera visualización de calidad QIIME2
        2. Crea gráficos de perfil de calidad con dokdo
        3. Opcionalmente filtra secuencias por q-score

        Args:
            output_dir: Directorio de salida
            filter_sequences: Si True, filtra las secuencias por calidad
            min_quality: Calidad mínima promedio para filtrado (por defecto: 4)

        Returns:
            dict: Diccionario con paths a los archivos generados
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print("🎯 Iniciando control de calidad...")

        results = {}

        # Paso 1: Visualización de calidad QIIME2
        print("\n📊 Generando visualización de calidad QIIME2...")
        viz_path = self.create_quality_visualization(output_dir)
        results['quality_viz'] = viz_path

        # Paso 2: Gráficos de calidad con dokdo
        if DOKDO_AVAILABLE:
            print("\n📈 Creando gráficos de perfil de calidad con dokdo...")
            plot_paths = self.plot_quality_profile_dokdo(output_path)
            results['quality_plots'] = plot_paths
        else:
            print("\n⚠️  dokdo no está disponible. Saltando generación de gráficos.")
            print("   Para instalar: pip install dokdo")
            results['quality_plots'] = None

        # Paso 3: Filtrado por q-score (opcional)
        if filter_sequences:
            print("\n🔬 Filtrando secuencias por q-score...")
            filtered_path, stats_path = self.filter_by_qscore(
                output_dir,
                min_quality=min_quality
            )
            results['filtered_seqs'] = filtered_path
            results['filter_stats'] = stats_path
        else:
            print("\n✓ Filtrado de secuencias omitido (use --filter-sequences para activar)")
            results['filtered_seqs'] = None
            results['filter_stats'] = None

        print("\n✅ Control de calidad completado exitosamente!")
        return results

    def create_quality_visualization(self, output_dir="results/quality"):
        """Crea la visualización de calidad de QIIME2"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print("  📊 Generando visualización de calidad...")
        self.quality_visualization = summarize(self.demux_seqs)

        # Guardar visualización
        viz_path = output_path / "demux_quality.qzv"
        self.quality_visualization.visualization.save(str(viz_path))
        self.quality_viz_path = str(viz_path)

        print(f"  ✓ Visualización guardada: {viz_path}")
        return viz_path

    def plot_quality_profile_dokdo(self, output_path):
        """
        Genera gráficos de perfil de calidad usando dokdo

        Args:
            output_path: Path object para el directorio de salida

        Returns:
            dict: Paths a los archivos PNG generados
        """
        if not DOKDO_AVAILABLE:
            print("  ⚠️  dokdo no disponible")
            return None

        # Asegurarse de que tenemos la visualización
        if self.quality_viz_path is None:
            self.create_quality_visualization(output_path)

        print("  📊 Generando gráficos de calidad con dokdo...")

        # Configurar estilo de seaborn
        sns.set_style("whitegrid")

        plot_paths = {}

        try:
            # Crear figura con 2 subplots (forward y reverse)
            fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(14, 6))

            # Plot forward reads
            print("    → Forward reads...")
            try:
                dokdo.read_quality_plot(
                    self.quality_viz_path,
                    strand='forward',
                    ax=ax1
                )
                ax1.set_title('Forward Read - Quality Profile',
                              fontsize=14, fontweight='bold')
                ax1.set_xlabel('Position (bp)', fontsize=12)
                ax1.set_ylabel('Quality Score', fontsize=12)
            except Exception as e:
                print(f"    ⚠️  Error plotting forward reads: {e}")
                ax1.text(0.5, 0.5, 'Forward reads\nnot available',
                         ha='center', va='center', fontsize=14,
                         transform=ax1.transAxes)

            # Plot reverse reads
            print("    → Reverse reads...")
            try:
                dokdo.read_quality_plot(
                    self.quality_viz_path,
                    strand='reverse',
                    ax=ax2
                )
                ax2.set_title('Reverse Read - Quality Profile',
                              fontsize=14, fontweight='bold')
                ax2.set_xlabel('Position (bp)', fontsize=12)
                ax2.set_ylabel('')  # Ocultar ylabel para el segundo gráfico
                ax2.set_yticklabels([])  # Ocultar etiquetas del eje Y
            except Exception as e:
                print(f"    ⚠️  Error plotting reverse reads: {e}")
                ax2.text(0.5, 0.5, 'Reverse reads\nnot available',
                         ha='center', va='center', fontsize=14,
                         transform=ax2.transAxes)

            # Ajustar diseño
            plt.tight_layout()

            # Guardar figura combinada
            combined_path = output_path / "quality_profile_combined.png"
            plt.savefig(combined_path, dpi=300, bbox_inches='tight')
            print(f"  ✓ Gráfico combinado guardado: {combined_path}")
            plot_paths['combined'] = str(combined_path)

            # También guardar como SVG
            combined_svg = output_path / "quality_profile_combined.svg"
            plt.savefig(combined_svg, bbox_inches='tight')
            print(f"  ✓ Gráfico SVG guardado: {combined_svg}")
            plot_paths['combined_svg'] = str(combined_svg)

            plt.close(fig)

            # Crear gráficos individuales
            # Forward solo
            fig_f, ax_f = plt.subplots(1, 1, figsize=(10, 6))
            try:
                dokdo.read_quality_plot(
                    self.quality_viz_path,
                    strand='forward',
                    ax=ax_f
                )
                ax_f.set_title('Forward Read - Quality Profile',
                               fontsize=14, fontweight='bold')
                ax_f.set_xlabel('Position (bp)', fontsize=12)
                ax_f.set_ylabel('Quality Score', fontsize=12)

                forward_path = output_path / "quality_profile_forward.png"
                plt.savefig(forward_path, dpi=300, bbox_inches='tight')
                print(f"  ✓ Gráfico forward guardado: {forward_path}")
                plot_paths['forward'] = str(forward_path)
            except Exception as e:
                print(f"  ⚠️  No se pudo crear gráfico forward: {e}")
            finally:
                plt.close(fig_f)

            # Reverse solo
            fig_r, ax_r = plt.subplots(1, 1, figsize=(10, 6))
            try:
                dokdo.read_quality_plot(
                    self.quality_viz_path,
                    strand='reverse',
                    ax=ax_r
                )
                ax_r.set_title('Reverse Read - Quality Profile',
                               fontsize=14, fontweight='bold')
                ax_r.set_xlabel('Position (bp)', fontsize=12)
                ax_r.set_ylabel('Quality Score', fontsize=12)

                reverse_path = output_path / "quality_profile_reverse.png"
                plt.savefig(reverse_path, dpi=300, bbox_inches='tight')
                print(f"  ✓ Gráfico reverse guardado: {reverse_path}")
                plot_paths['reverse'] = str(reverse_path)
            except Exception as e:
                print(f"  ⚠️  No se pudo crear gráfico reverse: {e}")
            finally:
                plt.close(fig_r)

        except Exception as e:
            print(f"  ❌ Error generando gráficos: {e}")
            import traceback
            print(f"  📋 Detalles: {traceback.format_exc()}")
            return None

        return plot_paths

    def filter_by_qscore(self, output_dir, min_quality=4):
        """
        Filtra secuencias por calidad usando q-score de QIIME2

        Args:
            output_dir: Directorio de salida
            min_quality: Calidad mínima promedio (por defecto: 4)

        Returns:
            tuple: (filtered_sequences_path, filter_stats_path)
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print(f"  🔬 Aplicando filtrado por q-score (min_quality={min_quality})...")

        try:
            # Aplicar filtrado q-score
            filter_result = q_score(
                demux=self.demux_seqs,
                min_quality=min_quality
            )

            # Guardar secuencias filtradas
            filtered_path = output_path / "filtered_seqs.qza"
            filter_result.filtered_sequences.save(str(filtered_path))
            print(f"  ✓ Secuencias filtradas guardadas: {filtered_path}")

            # Guardar estadísticas de filtrado
            stats_path = output_path / "filter_stats.qza"
            filter_result.filter_stats.save(str(stats_path))
            print(f"  ✓ Estadísticas de filtrado guardadas: {stats_path}")

            # Crear visualización de estadísticas usando metadata tabulate
            print("  📊 Generando visualización de estadísticas de filtrado...")
            try:
                from qiime2.plugins.metadata.visualizers import tabulate

                # Crear visualización directamente desde el artifact
                stats_viz = tabulate(filter_result.filter_stats)
                stats_viz_path = output_path / "filter_stats.qzv"
                stats_viz.visualization.save(str(stats_viz_path))
                print(f"  ✓ Visualización de stats guardada: {stats_viz_path}")
            except Exception as viz_error:
                print(f"  ⚠️  No se pudo crear visualización de stats: {viz_error}")
                print(f"  ℹ️  Puedes visualizar el archivo .qza directamente en QIIME2 View")

            # Generar gráficos de calidad de las secuencias FILTRADAS
            if DOKDO_AVAILABLE:
                print("\n  📈 Generando gráficos de calidad de secuencias FILTRADAS...")
                try:
                    # Crear visualización de calidad para las secuencias filtradas
                    print("    → Creando visualización de calidad QIIME2...")
                    filtered_quality_viz = summarize(filter_result.filtered_sequences)
                    filtered_viz_path = output_path / "filtered_quality.qzv"
                    filtered_quality_viz.visualization.save(str(filtered_viz_path))
                    print(f"    ✓ Visualización guardada: {filtered_viz_path}")

                    # Generar gráficos dokdo de las secuencias filtradas
                    filtered_plots = self._plot_filtered_quality(
                        str(filtered_viz_path),
                        output_path
                    )

                    if filtered_plots:
                        print(f"    ✓ Gráficos de secuencias filtradas generados:")
                        for plot_type, plot_path in filtered_plots.items():
                            print(f"      - {plot_type}: {plot_path}")

                except Exception as plot_error:
                    print(f"    ⚠️  Error generando gráficos de filtradas: {plot_error}")
            else:
                print("\n  ⚠️  dokdo no disponible, saltando gráficos de secuencias filtradas")

            # Mostrar resumen
            print("\n  📈 Resumen del filtrado:")
            print(f"     - Calidad mínima: {min_quality}")
            print(f"     - Secuencias filtradas: {filtered_path}")
            print(f"     - Estadísticas: {stats_path}")

            return str(filtered_path), str(stats_path)

        except Exception as e:
            print(f"  ❌ Error en filtrado por q-score: {e}")
            import traceback
            print(f"  📋 Detalles: {traceback.format_exc()}")
            raise

    def _plot_filtered_quality(self, filtered_viz_path, output_path):
        """
        Genera gráficos de perfil de calidad para secuencias filtradas usando dokdo

        Args:
            filtered_viz_path: Path al archivo .qzv de calidad de secuencias filtradas
            output_path: Path object para el directorio de salida

        Returns:
            dict: Paths a los archivos PNG generados
        """
        if not DOKDO_AVAILABLE:
            return None

        print("    📊 Generando gráficos con dokdo...")

        # Configurar estilo de seaborn
        sns.set_style("whitegrid")

        plot_paths = {}

        try:
            # Crear figura con 2 subplots (forward y reverse)
            fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(14, 6))

            # Plot forward reads
            print("      → Forward reads...")
            try:
                dokdo.read_quality_plot(
                    filtered_viz_path,
                    strand='forward',
                    ax=ax1
                )
                ax1.set_title('Filtered Forward Read - Quality Profile',
                              fontsize=14, fontweight='bold', color='darkgreen')
                ax1.set_xlabel('Position (bp)', fontsize=12)
                ax1.set_ylabel('Quality Score', fontsize=12)
            except Exception as e:
                print(f"      ⚠️  Error plotting forward reads: {e}")
                ax1.text(0.5, 0.5, 'Forward reads\nnot available',
                         ha='center', va='center', fontsize=14,
                         transform=ax1.transAxes)

            # Plot reverse reads
            print("      → Reverse reads...")
            try:
                dokdo.read_quality_plot(
                    filtered_viz_path,
                    strand='reverse',
                    ax=ax2
                )
                ax2.set_title('Filtered Reverse Read - Quality Profile',
                              fontsize=14, fontweight='bold', color='darkgreen')
                ax2.set_xlabel('Position (bp)', fontsize=12)
                ax2.set_ylabel('')  # Ocultar ylabel para el segundo gráfico
                ax2.set_yticklabels([])  # Ocultar etiquetas del eje Y
            except Exception as e:
                print(f"      ⚠️  Error plotting reverse reads: {e}")
                ax2.text(0.5, 0.5, 'Reverse reads\nnot available',
                         ha='center', va='center', fontsize=14,
                         transform=ax2.transAxes)

            # Ajustar diseño
            plt.tight_layout()

            # Guardar figura combinada con prefijo "filtered_"
            combined_path = output_path / "filtered_quality_profile_combined.png"
            plt.savefig(combined_path, dpi=300, bbox_inches='tight')
            plot_paths['combined'] = str(combined_path)

            # También guardar como SVG
            combined_svg = output_path / "filtered_quality_profile_combined.svg"
            plt.savefig(combined_svg, bbox_inches='tight')
            plot_paths['combined_svg'] = str(combined_svg)

            plt.close(fig)

            # Crear gráficos individuales
            # Forward solo
            fig_f, ax_f = plt.subplots(1, 1, figsize=(10, 6))
            try:
                dokdo.read_quality_plot(
                    filtered_viz_path,
                    strand='forward',
                    ax=ax_f
                )
                ax_f.set_title('Filtered Forward Read - Quality Profile',
                               fontsize=14, fontweight='bold', color='darkgreen')
                ax_f.set_xlabel('Position (bp)', fontsize=12)
                ax_f.set_ylabel('Quality Score', fontsize=12)

                forward_path = output_path / "filtered_quality_profile_forward.png"
                plt.savefig(forward_path, dpi=300, bbox_inches='tight')
                plot_paths['forward'] = str(forward_path)
            except Exception as e:
                print(f"      ⚠️  No se pudo crear gráfico forward: {e}")
            finally:
                plt.close(fig_f)

            # Reverse solo
            fig_r, ax_r = plt.subplots(1, 1, figsize=(10, 6))
            try:
                dokdo.read_quality_plot(
                    filtered_viz_path,
                    strand='reverse',
                    ax=ax_r
                )
                ax_r.set_title('Filtered Reverse Read - Quality Profile',
                               fontsize=14, fontweight='bold', color='darkgreen')
                ax_r.set_xlabel('Position (bp)', fontsize=12)
                ax_r.set_ylabel('Quality Score', fontsize=12)

                reverse_path = output_path / "filtered_quality_profile_reverse.png"
                plt.savefig(reverse_path, dpi=300, bbox_inches='tight')
                plot_paths['reverse'] = str(reverse_path)
            except Exception as e:
                print(f"      ⚠️  No se pudo crear gráfico reverse: {e}")
            finally:
                plt.close(fig_r)

        except Exception as e:
            print(f"    ❌ Error generando gráficos: {e}")
            import traceback
            print(f"    📋 Detalles: {traceback.format_exc()}")
            return None

        return plot_paths

    def get_visualization(self):
        """Retorna la visualización de calidad"""
        return self.quality_visualization

    def get_visualization_path(self):
        """Retorna el path a la visualización"""
        return self.quality_viz_path