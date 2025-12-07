"""
Módulo simplificado para control de calidad - Solo gráficos de calidad
"""
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from qiime2 import Artifact
from qiime2.plugins.demux.visualizers import summarize
import tempfile
import os


class QualityControl:
    """Clase simplificada para control de calidad"""

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
        self.quality_data = {}

    def run_quality_control(self, output_dir="results/quality_control"):
        """
        Ejecuta el control de calidad simplificado:
        1. Genera visualización de calidad QIIME2
        2. Crea gráficos de perfil de calidad

        Args:
            output_dir: Directorio de salida
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print("🎯 Iniciando control de calidad...")

        # Paso 1: Visualización de calidad QIIME2
        print("📊 Generando visualización de calidad QIIME2...")
        viz_path = self.create_quality_visualization(output_dir)

        # Paso 2: Extraer datos de calidad
        print("🔍 Extrayendo datos de calidad...")
        self.extract_quality_data()

        # Paso 3: Gráficos de calidad
        print("📈 Creando gráficos de perfil de calidad...")
        plot_path = self.plot_quality_profile(output_path / "quality_profile.png")

        print("✅ Control de calidad completado exitosamente!")

        return {
            'quality_viz': viz_path,
            'quality_plot': plot_path,
            'quality_data': self.quality_data
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

    def extract_quality_data(self):
        """Extrae datos de calidad de la visualización QIIME2"""
        if self.quality_visualization is None:
            self.create_quality_visualization()

        viz = self.quality_visualization.visualization

        with tempfile.TemporaryDirectory() as tmpdir:
            # Exportar datos
            viz.export_data(tmpdir)

            # Buscar archivos de calidad
            quality_files = {
                'forward': os.path.join(tmpdir, 'forward-seven-number-summaries.tsv'),
                'reverse': os.path.join(tmpdir, 'reverse-seven-number-summaries.tsv')
            }

            for strand, filepath in quality_files.items():
                if os.path.exists(filepath):
                    print(f"📄 Leyendo datos de calidad para {strand}...")
                    try:
                        # Leer el archivo TSV
                        df = pd.read_csv(filepath, sep='\t')
                        self.quality_data[strand] = df
                        print(f"✅ Datos de {strand} extraídos: {len(df)} posiciones")
                    except Exception as e:
                        print(f"⚠️  Error leyendo {strand}: {e}")
                else:
                    print(f"ℹ️  No se encontraron datos para {strand}")

    def plot_quality_profile(self, output_file="quality_profile.png"):
        """Genera gráficos de perfil de calidad"""
        if not self.quality_data:
            print("⚠️  No hay datos de calidad para graficar")
            return None

        print(f"📊 Generando gráfico de calidad...")

        # Determinar strands disponibles
        strands = list(self.quality_data.keys())
        num_strands = len(strands)

        if num_strands == 0:
            print("⚠️  No hay datos de calidad disponibles")
            return None

        # Configurar figura
        if num_strands == 2:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        else:
            fig, axes = plt.subplots(1, 1, figsize=(10, 6))
            axes = [axes]  # Convertir a lista para iterar

        colors = ['#1f77b4', '#ff7f0e']  # Colores para los gráficos

        for idx, strand in enumerate(strands):
            if idx < len(axes):
                ax = axes[idx]
                df = self.quality_data[strand]

                # Asegurarse de que tenemos las columnas necesarias
                required_cols = ['position', 'mean', 'median', 'min', 'max', 'q1', 'q3']
                available_cols = [col for col in required_cols if col in df.columns]

                if len(available_cols) < 2:
                    print(f"⚠️  Datos insuficientes para {strand}")
                    ax.text(0.5, 0.5, f'No hay datos\nde calidad\npara {strand}',
                            ha='center', va='center', fontsize=14,
                            transform=ax.transAxes)
                    ax.set_title(f'{strand.capitalize()} Read', fontsize=12)
                    continue

                # Extraer posiciones
                positions = df['position'].values if 'position' in df.columns else np.arange(len(df))

                # Graficar mediana
                if 'median' in df.columns:
                    ax.plot(positions, df['median'], color=colors[idx % 2],
                            linewidth=2, label='Mediana', zorder=5)

                # Graficar media
                if 'mean' in df.columns:
                    ax.plot(positions, df['mean'], color='#2ca02c',
                            linewidth=1.5, linestyle='--', label='Media', zorder=4)

                # Área entre percentiles (q1-q3)
                if all(col in df.columns for col in ['q1', 'q3']):
                    ax.fill_between(positions, df['q1'], df['q3'],
                                    color=colors[idx % 2], alpha=0.3, label='Q1-Q3')

                # Líneas de mínimo y máximo
                if all(col in df.columns for col in ['min', 'max']):
                    ax.fill_between(positions, df['min'], df['max'],
                                    color=colors[idx % 2], alpha=0.1, label='Min-Max')

                # Línea de referencia de calidad 20 y 30
                ax.axhline(y=20, color='red', linestyle=':', alpha=0.5, linewidth=1)
                ax.axhline(y=30, color='green', linestyle=':', alpha=0.5, linewidth=1)

                # Configurar ejes y etiquetas
                ax.set_xlabel('Posición en el read', fontsize=11)
                ax.set_ylabel('Score de Calidad', fontsize=11)
                ax.set_title(f'{strand.capitalize()} Read - Perfil de Calidad',
                             fontsize=12, fontweight='bold')

                # Añadir grid
                ax.grid(True, alpha=0.3, linestyle='--')

                # Añadir leyenda
                if idx == 0:  # Solo en el primer gráfico para no duplicar
                    ax.legend(loc='upper right', fontsize=9)

                # Ajustar límites del eje Y
                y_min = 0
                if 'min' in df.columns:
                    y_min = min(y_min, df['min'].min() * 0.9)
                ax.set_ylim(bottom=y_min, top=45)

        plt.tight_layout()

        # Guardar figura
        try:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✅ Gráfico de calidad guardado: {output_file}")

            # También guardar una versión en formato vectorial
            output_svg = output_file.with_suffix('.svg')
            plt.savefig(output_svg, bbox_inches='tight')
            print(f"✅ Gráfico de calidad (SVG) guardado: {output_svg}")

        except Exception as e:
            print(f"❌ Error guardando gráfico: {e}")
            return None
        finally:
            plt.close(fig)

        return output_file

    def get_visualization(self):
        """Retorna la visualización de calidad"""
        return self.quality_visualization

    def get_quality_stats(self):
        """Retorna estadísticas de calidad"""
        stats = {}
        for strand, df in self.quality_data.items():
            if not df.empty:
                stats[strand] = {
                    'num_positions': len(df),
                    'mean_quality': df['mean'].mean() if 'mean' in df.columns else None,
                    'median_quality': df['median'].median() if 'median' in df.columns else None,
                    'min_quality': df['min'].min() if 'min' in df.columns else None,
                    'max_quality': df['max'].max() if 'max' in df.columns else None,
                    'positions_below_20': len(df[df['mean'] < 20]) if 'mean' in df.columns else None,
                    'positions_below_30': len(df[df['mean'] < 30]) if 'mean' in df.columns else None,
                }
        return stats