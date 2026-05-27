# -*- coding: utf-8 -*-
"""
Created on Wed May 20 23:02:06 2026

@author: user
"""

from pathlib import Path
import pandas as pd
import networkx as nx


def normalize_node_id(node_id: str, source: str = "kegg") -> str:
    node_id = str(node_id).strip()

    if source.lower() == "kegg":
        if node_id.startswith("KEGG:"):
            return node_id
        if node_id.startswith("hsa:"):
            return "KEGG:" + node_id.replace("hsa:", "")
        return "KEGG:" + node_id

    if source.lower() == "reactome":
        if node_id.startswith("REACTOME:"):
            return node_id
        return "REACTOME:" + node_id

    return node_id


def load_graph_from_csv(results_dir: str) -> nx.DiGraph:
    results_dir = Path(results_dir)

    nodes_path = results_dir / "nodes.csv"
    edges_path = results_dir / "edges.csv"

    if not nodes_path.exists():
        raise FileNotFoundError(f"nodes.csv not found: {nodes_path}")

    if not edges_path.exists():
        raise FileNotFoundError(f"edges.csv not found: {edges_path}")

    nodes_df = pd.read_csv(nodes_path)
    edges_df = pd.read_csv(edges_path)

    graph = nx.DiGraph()

    for _, row in nodes_df.iterrows():
        graph.add_node(
            row["node_id"],
            source=row.get("source", ""),
            type=row.get("type", ""),
            label=row.get("label", ""),
            raw_id=row.get("raw_id", ""),
        )

    for _, row in edges_df.iterrows():
        graph.add_edge(
            row["source_node"],
            row["target_node"],
            relation_type=row.get("relation_type", ""),
            source=row.get("source", ""),
        )

    return graph


def check_path(results_dir: str, start_node: str, end_node: str, source: str = "kegg", cutoff: int = 10):
    graph = load_graph_from_csv(results_dir)

    start_node = normalize_node_id(start_node, source)
    end_node = normalize_node_id(end_node, source)

    print(f"Start node: {start_node}")
    print(f"End node  : {end_node}")

    if start_node not in graph:
        print(f"시작 노드가 nodes.csv에 없습니다: {start_node}")
        return

    if end_node not in graph:
        print(f"도착 노드가 nodes.csv에 없습니다: {end_node}")
        return

    print(f"Total nodes: {graph.number_of_nodes()}")
    print(f"Total edges: {graph.number_of_edges()}")

    if nx.has_path(graph, start_node, end_node):
        paths = list(nx.all_simple_paths(
            graph,
            source=start_node,
            target=end_node,
            cutoff=cutoff
        ))

        print(f"\n방향성 경로가 있습니다. 총 {len(paths)}개 발견.")

        for i, path in enumerate(paths, start=1):
            path_labels = []
            for node in path:
                label = graph.nodes[node].get("label", "")
                path_labels.append(f"{node}({label})")

            print(f"[Path {i}] " + " -> ".join(path_labels))

    else:
        print("\n방향성 경로는 없습니다.")

        reverse_exists = nx.has_path(graph, end_node, start_node)
        if reverse_exists:
            print(f"단, 반대 방향 경로는 존재합니다: {end_node} -> {start_node}")

        undirected_graph = graph.to_undirected()
        if nx.has_path(undirected_graph, start_node, end_node):
            print("무방향 연결성은 존재합니다.")
            shortest = nx.shortest_path(undirected_graph, start_node, end_node)

            print("\n무방향 최단 연결:")
            print(" -> ".join(
                f"{node}({graph.nodes[node].get('label', '')})"
                for node in shortest
            ))
        else:
            print("무방향으로도 두 노드는 연결되어 있지 않습니다.")


if __name__ == "__main__":
    check_path(
        results_dir=r"results\kegg_hsa04066_20260520_222849",
        start_node="3091",
        end_node="7422",
        source="kegg",
        cutoff=10,
    )