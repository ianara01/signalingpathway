# -*- coding: utf-8 -*-
"""
Created on Mon May 25 00:57:52 2026

@author: user
"""

from pathlib import Path
from datetime import datetime
import csv
import json

import networkx as nx

from SignalPathway.s2pathway import SignalingPathway

    
def save_graph_tables(graph: nx.DiGraph, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    nodes_path = output_dir / "nodes.csv"
    graph_json_path = output_dir / "graph_node_link.json"
    edges_path = output_dir / "edges.csv"
    graph_json_path = output_dir / "graph_edge_link.json"

    with nodes_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "node_id",
            "source",
            "type",
            "label",
            "raw_id",
            "compound_id",
            "compound_name",
            "kgml_entry_id",
            "kgml_name",
            "group_components",
            "group_component_count",
        ])

        for node_id, attrs in graph.nodes(data=True):
            writer.writerow([
                node_id,
                attrs.get("source", ""),
                attrs.get("type", ""),
                attrs.get("label", ""),
                attrs.get("raw_id", ""),
                attrs.get("compound_id", ""),
                attrs.get("compound_name", ""),
                attrs.get("kgml_entry_id", ""),
                attrs.get("kgml_name", ""),
                attrs.get("group_components", ""),
                attrs.get("group_component_count", ""),
            ])

    with edges_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "source_node",
            "target_node",
            "relation_type",
            "source",
            "kgml_relation_type",
        ])

        for u, v, attrs in graph.edges(data=True):
            writer.writerow([
                u,
                v,
                attrs.get("relation_type", ""),
                attrs.get("source", ""),
                attrs.get("kgml_relation_type", ""),
            ])

    with graph_json_path.open("w", encoding="utf-8") as f:
        json.dump(nx.node_link_data(graph), f, ensure_ascii=False, indent=2)


def save_node_filter(graph: nx.DiGraph, output_path: Path, allowed_types: set[str]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "source", "type", "label", "raw_id"])

        for node_id, attrs in graph.nodes(data=True):
            if attrs.get("type") in allowed_types:
                writer.writerow([
                    node_id,
                    attrs.get("source", ""),
                    attrs.get("type", ""),
                    attrs.get("label", ""),
                    attrs.get("raw_id", ""),
                ])


def save_edge_filter(graph: nx.DiGraph, output_path: Path, relation_keywords: list[str]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["source_node", "target_node", "relation_type", "source"])

        for u, v, attrs in graph.edges(data=True):
            relation_type = str(attrs.get("relation_type", ""))
            if any(keyword.lower() in relation_type.lower() for keyword in relation_keywords):
                writer.writerow([
                    u,
                    v,
                    relation_type,
                    attrs.get("source", ""),
                ])


def find_nodes_by_keywords(graph: nx.DiGraph, keywords: list[str]) -> dict[str, list[str]]:
    results = {}

    for keyword in keywords:
        query = keyword.upper()
        matched = []

        for node_id, attrs in graph.nodes(data=True):
            haystack = " ".join([
                str(node_id),
                str(attrs.get("label", "")),
                str(attrs.get("raw_id", "")),
                str(attrs.get("compound_id", "")),
                str(attrs.get("compound_name", "")),
                str(attrs.get("kgml_name", "")),
                str(attrs.get("group_components", "")),
            ]).upper()
            
            if query in haystack:
                matched.append(node_id)

        results[keyword] = matched

    return results


    def save_keyword_hits(graph: nx.DiGraph, keyword_hits: dict[str, list[str]], output_path: Path) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
        with output_path.open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "keyword",
                "node_id",
                "label",
                "type",
                "source",
                "raw_id",
                "compound_id",
                "compound_name",
                "kgml_entry_id",
                "kgml_name",
                "group_components",
                "group_component_count",
            ])
    
            for keyword, node_ids in keyword_hits.items():
                for node_id in node_ids:
                    attrs = graph.nodes[node_id]
                    writer.writerow([
                        keyword,
                        node_id,
                        attrs.get("label", ""),
                        attrs.get("type", ""),
                        attrs.get("source", ""),
                        attrs.get("raw_id", ""),
                        attrs.get("compound_id", ""),
                        attrs.get("compound_name", ""),
                        attrs.get("kgml_entry_id", ""),
                        attrs.get("kgml_name", ""),
                        attrs.get("group_components", ""),
                        attrs.get("group_component_count", ""),
                    ])
                    

def save_paths(
    graph: nx.DiGraph,
    start_node: str,
    end_node: str,
    output_path: Path,
    cutoff: int = 10,
    undirected: bool = False,
) -> list[list[str]]:
    output_path.parent.mkdir(parents=True, exist_ok=True)

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

def save_compounds(graph: nx.DiGraph, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "node_id",
            "compound_id",
            "compound_name",
            "label",
            "type",
            "raw_id",
        ])

        for node_id, attrs in graph.nodes(data=True):
            if attrs.get("type") != "compound":
                continue

            compound_id = (
                attrs.get("compound_id")
                or attrs.get("raw_id")
                or str(node_id).replace("KEGG:", "")
            )

            writer.writerow([
                node_id,
                compound_id,
                attrs.get("compound_name", ""),
                attrs.get("label", ""),
                attrs.get("type", ""),
                attrs.get("raw_id", ""),
            ])

def save_groups(graph: nx.DiGraph, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "node_id",
            "label",
            "raw_id",
            "kgml_entry_id",
            "kgml_name",
            "group_components",
            "group_component_count",
        ])

        for node_id, attrs in graph.nodes(data=True):
            if attrs.get("type") != "group":
                continue

            writer.writerow([
                node_id,
                attrs.get("label", ""),
                attrs.get("raw_id", ""),
                attrs.get("kgml_entry_id", ""),
                attrs.get("kgml_name", ""),
                attrs.get("group_components", ""),
                attrs.get("group_component_count", ""),
            ])

def save_keyword_hits(graph: nx.DiGraph, keyword_hits: dict[str, list[str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "keyword",
            "node_id",
            "label",
            "type",
            "source",
            "raw_id",
            "compound_id",
            "compound_name",
            "kgml_entry_id",
            "kgml_name",
            "group_components",
            "group_component_count",
        ])

        for keyword, node_ids in keyword_hits.items():
            for node_id in node_ids:
                attrs = graph.nodes[node_id]
                writer.writerow([
                    keyword,
                    node_id,
                    attrs.get("label", ""),
                    attrs.get("type", ""),
                    attrs.get("source", ""),
                    attrs.get("raw_id", ""),
                    attrs.get("compound_id", ""),
                    attrs.get("compound_name", ""),
                    attrs.get("kgml_entry_id", ""),
                    attrs.get("kgml_name", ""),
                    attrs.get("group_components", ""),
                    attrs.get("group_component_count", ""),
                ])
                
def run_kegg_hif1_workflow() -> Path:
    """
    KEGG HIF-1 workflow.

    KEGG에서는 Reactome처럼 Pathway/Entity/Reaction을 따로 조합하기보다,
    HIF-1 signaling pathway인 hsa04066 KGML을 가져와 gene/compound/relation을 분석한다.
    """
    pathway_id = "hsa04066"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    output_dir = Path("results") / f"kegg_hif1_{pathway_id}_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    pathway = SignalingPathway(source="kegg")
    result = pathway.fetch_pathway(pathway_id)

    graph = pathway.graph

    save_graph_tables(graph, output_dir)

    # KEGG HIF-1 핵심 node 검색
    keyword_hits = find_nodes_by_keywords(
        graph,
        keywords=[
            "HIF1A",
            "HIF-1",
            "VEGF",
            "VEGFA",
            "MTOR",
            "mTOR",
            "ARNT",
            "VHL",
            "EGLN",
            "PHD",
            "PDK1",
            "SLC2A1",
        ],
    )
    save_keyword_hits(graph, keyword_hits, output_dir / "kegg_hif1_keyword_hits.csv")

    # KEGG compound node 추출
    save_compounds(
        graph,
        output_dir / "kegg_compounds.csv",
    )

    # KEGG activation/inhibition/indirect effect edge 분리
    save_edge_filter(
        graph,
        output_dir / "kegg_activation_edges.csv",
        relation_keywords=["activation"],
    )

    save_edge_filter(
        graph,
        output_dir / "kegg_inhibition_edges.csv",
        relation_keywords=["inhibition"],
    )

    save_edge_filter(
        graph,
        output_dir / "kegg_indirect_effect_edges.csv",
        relation_keywords=["indirect effect"],
    )
    save_groups(
        graph,
        output_dir / "kegg_groups.csv",
    )

    # 대표 경로 확인: HIF1A -> VEGFA
    # KEGG node id는 보통 KEGG:3091 = HIF1A, KEGG:7422 = VEGFA
    directed_paths = save_paths(
        graph,
        start_node="KEGG:3091",
        end_node="KEGG:7422",
        output_path=output_dir / "paths_directed_HIF1A_to_VEGFA.csv",
        cutoff=10,
        undirected=False,
    )

    undirected_paths = save_paths(
        graph,
        start_node="KEGG:3091",
        end_node="KEGG:7422",
        output_path=output_dir / "paths_undirected_HIF1A_to_VEGFA.csv",
        cutoff=10,
        undirected=True,
    )

    directed_exists = (
        "KEGG:3091" in graph
        and "KEGG:7422" in graph
        and nx.has_path(graph, "KEGG:3091", "KEGG:7422")
    )
    
    reverse_exists = (
        "KEGG:3091" in graph
        and "KEGG:7422" in graph
        and nx.has_path(graph, "KEGG:7422", "KEGG:3091")
    )
    
    undirected_exists = (
        "KEGG:3091" in graph
        and "KEGG:7422" in graph
        and nx.has_path(graph.to_undirected(), "KEGG:3091", "KEGG:7422")
    )
    
    summary = {
        "workflow": "kegg_hif1",
        "pathway_id": pathway_id,
        "result": result,
        "node_count": graph.number_of_nodes(),
        "edge_count": graph.number_of_edges(),
        "HIF1A_node_exists": "KEGG:3091" in graph,
        "VEGFA_node_exists": "KEGG:7422" in graph,
        "directed_HIF1A_to_VEGFA_exists": directed_exists,
        "reverse_VEGFA_to_HIF1A_exists": reverse_exists,
        "undirected_HIF1A_to_VEGFA_exists": undirected_exists,
        "directed_HIF1A_to_VEGFA_path_count": len(directed_paths),
        "undirected_HIF1A_to_VEGFA_path_count": len(undirected_paths),
        "output_dir": str(output_dir),
    }

    with (output_dir / "workflow_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    pathway.visualize_nodesedges(
        with_edge_labels=False,
        save_path=output_dir / "kegg_hif1_graph.png",
        show=False,
    )

    print(f"KEGG HIF-1 workflow saved: {output_dir}")
    return output_dir


if __name__ == "__main__":
    run_kegg_hif1_workflow()