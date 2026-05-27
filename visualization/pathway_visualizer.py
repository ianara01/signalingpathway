"""Pathway visualization utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx


DEFAULT_COLOR_MAP = {
    # KEGG
    "gene": "lightblue",
    "compound": "orange",
    "map": "lightgray",
    "group": "plum",
    # Reactome
    "Protein": "lightblue",
    "EntityWithAccessionedSequence": "lightblue",
    "ChemicalCompound": "orange",
    "SimpleEntity": "orange",
    "Complex": "plum",
    "Reaction": "lightgreen",
    "ReactionLikeEvent": "lightgreen",
    "Pathway": "lightgray",
}


def visualize_graph(
    graph,
    title=None,
    with_edge_labels=True,
    show=True,
    save_path=None,
):
    if not graph.nodes:
        print("시각화할 데이터가 없습니다.")
        return

    plt.figure(figsize=(14, 10))

    pos = nx.spring_layout(graph, seed=42, k=0.5)

    color_map = {
        "gene": "lightblue",
        "Protein": "lightblue",
        "compound": "orange",
        "ChemicalCompound": "orange",
        "Reaction": "lightgreen",
        "Complex": "plum",
        "map": "lightgray",
    }

    node_colors = [
        color_map.get(graph.nodes[node].get("type"), "gray")
        for node in graph.nodes
    ]

    labels = {
        node: graph.nodes[node].get("label", node)
        for node in graph.nodes
    }

    nx.draw(
        graph,
        pos,
        labels=labels,
        with_labels=True,
        node_color=node_colors,
        node_size=1400,
        font_size=8,
        edge_color="gray",
        arrows=True,
        arrowsize=12,
    )

    if with_edge_labels:
        edge_labels = nx.get_edge_attributes(graph, "relation_type")
        if edge_labels:
            nx.draw_networkx_edge_labels(
                graph,
                pos,
                edge_labels=edge_labels,
                font_size=7,
            )

    if title is None:
        source = graph.graph.get("source", "unknown")
        pathway_id = graph.graph.get("pathway_id", "unknown")
        title = f"{source.upper()} Pathway: {pathway_id}"

    plt.title(title)

    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)

        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Graph image saved: {save_path}")

    if show:
        plt.show()
    else:
        plt.close()


class PathwayVisualizer:
    """Small wrapper retained for compatibility with the previous code style."""

    def __init__(self, graph: nx.Graph):
        self.graph = graph

    def visualize_pathway(self, **kwargs: Any) -> None:
        visualize_graph(self.graph, **kwargs)
