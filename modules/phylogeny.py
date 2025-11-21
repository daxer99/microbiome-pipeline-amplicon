# modules/phylogeny.py
import os
from qiime2 import Artifact
from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree


def make_phylogeny(rep_seqs, output_folder):
    """Generar árbol filogenético a partir de secuencias representativas

    Args:
        rep_seqs: Ruta al artefacto QIIME2 de secuencias representativas o el artefacto mismo
        output_folder: Directorio donde guardar los resultados

    Returns:
        tuple: Rutas a los árboles sin raíz y con raíz
    """
    os.makedirs(output_folder, exist_ok=True)

    # Cargar artefacto si se pasa como string
    if isinstance(rep_seqs, str):
        rep_seqs = Artifact.load(rep_seqs)

    # Generar filogenia
    phylo = align_to_tree_mafft_fasttree(sequences=rep_seqs, n_threads='auto')

    # Guardar árboles
    unrooted_tree = phylo.tree
    rooted_tree = phylo.rooted_tree

    unrooted_tree.save(f"{output_folder}/unrooted_tree.qza")
    rooted_tree.save(f"{output_folder}/rooted_tree.qza")

    return f"{output_folder}/unrooted_tree.qza", f"{output_folder}/rooted_tree.qza"