# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:11:50 2026

@author: user
"""

# workflows/generic_kegg_target_workflow.py

from pathlib import Path
from datetime import datetime
import csv
import json

import networkx as nx

from SignalPathway.s2pathway import SignalingPathway


def save_graph_tables(graph: nx.DiGraph, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    with (output_dir / "nodes.csv").open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "source", "type", "label", "raw_id"])

        for node_id, attrs in graph.nodes(data=True):
            writer.writerow([
                node_id,
                attrs.get("source", ""),
                attrs.get("type", ""),
                attrs.get("label", ""),
                attrs.get("raw_id", ""),
            ])

    with (output_dir / "edges.csv").open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["source_node", "target_node", "relation_type", "source"])

        for u, v, attrs in graph.edges(data=True):
            writer.writerow([
                u,
                v,
                attrs.get("relation_type", attrs.get("relation", attrs.get("rel", ""))),
                attrs.get("source", ""),
            ])

    with (output_dir / "graph_node_link.json").open("w", encoding="utf-8") as f:
        json.dump(nx.node_link_data(graph), f, ensure_ascii=False, indent=2)


def find_nodes_by_keywords(graph: nx.DiGraph, keywords: list[str]) -> dict[str, list[str]]:
    results = {}

    for keyword in keywords:
        query = keyword.upper()
        hits = []

        for node_id, attrs in graph.nodes(data=True):
            haystack = " ".join([
                node_id,
                str(attrs.get("label", "")),
                str(attrs.get("raw_id", "")),
            ]).upper()

            if query in haystack:
                hits.append(node_id)

        results[keyword] = hits

    return results


def save_keyword_hits(graph: nx.DiGraph, hits: dict[str, list[str]], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["keyword", "node_id", "label", "type", "source", "raw_id"])

        for keyword, node_ids in hits.items():
            for node_id in node_ids:
                attrs = graph.nodes[node_id]
                writer.writerow([
                    keyword,
                    node_id,
                    attrs.get("label", ""),
                    attrs.get("type", ""),
                    attrs.get("source", ""),
                    attrs.get("raw_id", ""),
                ])


def save_compounds(graph: nx.DiGraph, output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "label", "type", "raw_id"])

        for node_id, attrs in graph.nodes(data=True):
            if attrs.get("type") == "compound":
                writer.writerow([
                    node_id,
                    attrs.get("label", ""),
                    attrs.get("type", ""),
                    attrs.get("raw_id", ""),
                ])


def save_relation_edges(graph: nx.DiGraph, output_dir: Path) -> None:
    groups = {
        "activation_edges.csv": ["activation"],
        "inhibition_edges.csv": ["inhibition"],
        "indirect_effect_edges.csv": ["indirect effect"],
        "expression_edges.csv": ["expression"],
        "phosphorylation_edges.csv": ["phosphorylation"],
    }

    for filename, keywords in groups.items():
        with (output_dir / filename).open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["source_node", "target_node", "relation_type", "source"])

            for u, v, attrs in graph.edges(data=True):
                relation_type = str(attrs.get("relation_type", ""))

                if any(k.lower() in relation_type.lower() for k in keywords):
                    writer.writerow([
                        u,
                        v,
                        relation_type,
                        attrs.get("source", ""),
                    ])


def save_paths(
    graph: nx.DiGraph,
    start_node: str,
    end_node: str,
    output_path: Path,
    cutoff: int = 10,
    undirected: bool = False,
) -> list[list[str]]:
    search_graph = graph.to_undirected() if undirected else graph

    if start_node not in search_graph or end_node not in search_graph:
        paths = []
    elif nx.has_path(search_graph, start_node, end_node):
        paths = list(nx.all_simple_paths(search_graph, start_node, end_node, cutoff=cutoff))
    else:
        paths = []

    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["path_index", "step_index", "node_id", "label", "type"])

        for path_index, path in enumerate(paths, start=1):
            for step_index, node_id in enumerate(path, start=1):
                attrs = graph.nodes[node_id]
                writer.writerow([
                    path_index,
                    step_index,
                    node_id,
                    attrs.get("label", ""),
                    attrs.get("type", ""),
                ])

    return paths


def run_generic_kegg_target_workflow(target_key: str, target_config: dict) -> Path:
    kegg_config = target_config["kegg"]

    pathway_id = kegg_config["pathway_id"]
    pathway_name = kegg_config.get("pathway_name", pathway_id)
    keywords = kegg_config.get("keywords", [])
    representative_paths = kegg_config.get("representative_paths", [])

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path("results") / f"kegg_{target_key}_{pathway_id}_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    pathway = SignalingPathway(source="kegg")
    result = pathway.fetch_pathway(pathway_id)
    graph = pathway.graph

    save_graph_tables(graph, output_dir)

    hits = find_nodes_by_keywords(graph, keywords)
    save_keyword_hits(graph, hits, output_dir / f"kegg_{target_key}_keyword_hits.csv")

    save_compounds(graph, output_dir / "kegg_compounds.csv")
    save_relation_edges(graph, output_dir)

    path_summaries = []

    for path_cfg in representative_paths:
        name = path_cfg["name"]
        start = path_cfg["start"]
        end = path_cfg["end"]

        directed_paths = save_paths(
            graph,
            start_node=start,
            end_node=end,
            output_path=output_dir / f"paths_directed_{name}.csv",
            cutoff=10,
            undirected=False,
        )

        undirected_paths = save_paths(
            graph,
            start_node=start,
            end_node=end,
            output_path=output_dir / f"paths_undirected_{name}.csv",
            cutoff=10,
            undirected=True,
        )

        path_summaries.append({
            "name": name,
            "start": start,
            "end": end,
            "directed_path_count": len(directed_paths),
            "undirected_path_count": len(undirected_paths),
        })

    summary = {
        "workflow": "generic_kegg_target",
        "target_key": target_key,
        "pathway_id": pathway_id,
        "pathway_name": pathway_name,
        "node_count": graph.number_of_nodes(),
        "edge_count": graph.number_of_edges(),
        "path_summaries": path_summaries,
        "result": result,
        "output_dir": str(output_dir),
    }

    with (output_dir / "workflow_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    try:
        pathway.visualize_nodesedges(
            with_edge_labels=False,
            save_path=output_dir / f"kegg_{target_key}_graph.png",
            show=False,
        )
    except Exception as exc:
        print(f"graph image skipped: {exc}")

    print(f"KEGG {target_key} workflow saved: {output_dir}")
    return output_dir