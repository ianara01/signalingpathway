# -*- coding: utf-8 -*-
"""
Created on Mon May 25 00:58:58 2026

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
    edges_path = output_dir / "edges.csv"
    graph_json_path = output_dir / "graph_node_link.json"

    with nodes_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "node_id",
            "source",
            "type",
            "label",
            "raw_id",
            "schemaClass",
            "compartment",
            "species",
            "reference_identifier",
            "reference_database",
            "reference_display",
        ])

        for node_id, attrs in graph.nodes(data=True):
            writer.writerow([
                node_id,
                attrs.get("source", ""),
                attrs.get("type", ""),
                attrs.get("label", ""),
                attrs.get("raw_id", ""),
                attrs.get("schemaClass", ""),
                attrs.get("compartment", ""),
                attrs.get("species", ""),
                attrs.get("reference_identifier", ""),
                attrs.get("reference_database", ""),
                attrs.get("reference_display", ""),
            ])

    with edges_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "source_node",
            "target_node",
            "relation_type",
            "source",
        ])

        for u, v, attrs in graph.edges(data=True):
            writer.writerow([
                u,
                v,
                attrs.get("relation_type", attrs.get("relation", attrs.get("rel", ""))),
                attrs.get("source", ""),
            ])

    with graph_json_path.open("w", encoding="utf-8") as f:
        json.dump(nx.node_link_data(graph), f, ensure_ascii=False, indent=2)


def save_relation_edges(graph: nx.DiGraph, output_dir: Path) -> None:
    relation_groups = {
        "pathway_event_links.csv": {"has_event"},
        "reaction_io.csv": {"input", "output"},
        "catalyst_regulator_links.csv": {
            "catalyst",
            "regulator",
            "PositiveRegulation",
            "NegativeRegulation",
            "Requirement",
        },
        "entity_complex_links.csv": {"component_of", "component", "member", "candidate"},
        "entity_pathway_links.csv": {"located_in_pathway"},
        "entity_event_links.csv": {"participates_in_candidate_event"},
        "event_order_links.csv": {"precedingEvent"},
    }

    for filename, relation_types in relation_groups.items():
        output_path = output_dir / filename

        with output_path.open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "source_node",
                "source_label",
                "source_type",
                "target_node",
                "target_label",
                "target_type",
                "relation_type",
            ])

            for u, v, attrs in graph.edges(data=True):
                relation_type = attrs.get("relation_type", "")

                if relation_type not in relation_types:
                    continue

                writer.writerow([
                    u,
                    graph.nodes[u].get("label", ""),
                    graph.nodes[u].get("type", ""),
                    v,
                    graph.nodes[v].get("label", ""),
                    graph.nodes[v].get("type", ""),
                    relation_type,
                ])


def save_node_type_filter(
    graph: nx.DiGraph,
    output_path: Path,
    allowed_types: set[str],
) -> None:
    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "node_id",
            "label",
            "type",
            "source",
            "raw_id",
            "compartment",
            "species",
            "reference_identifier",
            "reference_database",
            "reference_display",
        ])

        for node_id, attrs in graph.nodes(data=True):
            if attrs.get("type") not in allowed_types:
                continue

            writer.writerow([
                node_id,
                attrs.get("label", ""),
                attrs.get("type", ""),
                attrs.get("source", ""),
                attrs.get("raw_id", ""),
                attrs.get("compartment", ""),
                attrs.get("species", ""),
                attrs.get("reference_identifier", ""),
                attrs.get("reference_database", ""),
                attrs.get("reference_display", ""),
            ])


def find_nodes_by_keywords(graph: nx.DiGraph, keywords: list[str]) -> dict[str, list[str]]:
    results = {}

    for keyword in keywords:
        query = keyword.upper()
        matched = []

        for node_id, attrs in graph.nodes(data=True):
            haystack = " ".join([
                node_id,
                str(attrs.get("label", "")),
                str(attrs.get("raw_id", "")),
                str(attrs.get("reference_identifier", "")),
                str(attrs.get("reference_display", "")),
            ]).upper()

            if query in haystack:
                matched.append(node_id)

        results[keyword] = matched

    return results


def save_keyword_hits(graph: nx.DiGraph, keyword_hits: dict[str, list[str]], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "keyword",
            "node_id",
            "label",
            "type",
            "raw_id",
            "compartment",
            "reference_identifier",
            "reference_display",
        ])

        for keyword, node_ids in keyword_hits.items():
            for node_id in node_ids:
                attrs = graph.nodes[node_id]
                writer.writerow([
                    keyword,
                    node_id,
                    attrs.get("label", ""),
                    attrs.get("type", ""),
                    attrs.get("raw_id", ""),
                    attrs.get("compartment", ""),
                    attrs.get("reference_identifier", ""),
                    attrs.get("reference_display", ""),
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
        writer.writerow(["path_index", "step_index", "node_id", "label", "type", "relation_to_next"])

        for path_index, path in enumerate(paths, start=1):
            for step_index, node_id in enumerate(path, start=1):
                relation_to_next = ""

                if step_index < len(path):
                    next_node = path[step_index]
                    if graph.has_edge(node_id, next_node):
                        relation_to_next = graph.edges[node_id, next_node].get("relation_type", "")
                    elif graph.has_edge(next_node, node_id):
                        relation_to_next = "reverse:" + graph.edges[next_node, node_id].get("relation_type", "")

                attrs = graph.nodes[node_id]
                writer.writerow([
                    path_index,
                    step_index,
                    node_id,
                    attrs.get("label", ""),
                    attrs.get("type", ""),
                    relation_to_next,
                ])

    return paths


def fetch_and_merge(
    master_graph: nx.DiGraph,
    reactome_id: str,
    *,
    build_reactome_flow: bool,
    recursive_events: bool,
    recursive_max_depth: int,
    step_output_dir: Path,
) -> dict:
    """
    SignalingPathway.fetch_pathway()는 매번 graph를 reset하므로,
    각 ID 실행 결과를 master_graph에 compose하여 통합한다.
    """
    pathway = SignalingPathway(source="reactome", timeout=(20, 180))

    result = pathway.fetch_pathway(
        reactome_id,
        build_reactome_flow=build_reactome_flow,
        recursive_events=recursive_events,
        recursive_max_depth=recursive_max_depth,
    )

    # 핵심 추가:
    # save_graph_tables() 전에 Reactome node 속성 보강
    enrich_reactome_node_attributes(pathway)

    step_output_dir.mkdir(parents=True, exist_ok=True)
    save_graph_tables(pathway.graph, step_output_dir)

    try:
        pathway.visualize_nodesedges(
            with_edge_labels=False,
            save_path=step_output_dir / "graph.png",
            show=False,
        )
    except Exception as exc:
        print(f"graph image skipped for {reactome_id}: {exc}")

    # enrichment 이후의 graph를 master_graph에 병합
    master_graph.update(pathway.graph)

    return result

def enrich_reactome_node_attributes(pathway):
    """
    Reactome graph의 node 중 속성이 부족한 node를 /query/{raw_id}로 다시 조회하여
    schemaClass, compartment, species, reference 정보를 보강한다.
    """
    if pathway.graph.graph.get("source") != "reactome":
        return

    for node_id, attrs in list(pathway.graph.nodes(data=True)):
        raw_id = attrs.get("raw_id")

        if not raw_id:
            continue

        raw_id = str(raw_id)

        # Reactome stable ID만 조회
        if not raw_id.startswith("R-"):
            continue

        # schemaClass만 있다고 충분한 것이 아니므로,
        # compartment/species/reference 계열 중 하나라도 비어 있으면 보강 시도
        needs_enrich = any([
            not attrs.get("schemaClass"),
            not attrs.get("compartment"),
            not attrs.get("species"),
            not attrs.get("reference_identifier"),
            not attrs.get("reference_database"),
            not attrs.get("reference_display"),
        ])

        if not needs_enrich:
            continue

        try:
            detail = pathway.reactome_provider.query_reactome_id(raw_id)
        except Exception:
            continue

        attrs["schemaClass"] = detail.get("schemaClass", attrs.get("schemaClass", ""))
        attrs["type"] = detail.get("schemaClass", attrs.get("type", ""))
        attrs["label"] = detail.get("displayName", attrs.get("label", raw_id))

        # compartment
        compartments = detail.get("compartment") or detail.get("compartments") or []
        if isinstance(compartments, dict):
            compartments = [compartments]

        for comp in compartments:
            if isinstance(comp, dict):
                attrs["compartment"] = comp.get("displayName", comp.get("name", ""))
                break

        # species
        species = detail.get("species")
        if isinstance(species, dict):
            attrs["species"] = species.get("displayName", species.get("name", ""))

        # referenceEntity
        ref = detail.get("referenceEntity")
        if isinstance(ref, dict):
            attrs["reference_identifier"] = ref.get("identifier", "")
            attrs["reference_database"] = ref.get("databaseName", "")
            attrs["reference_display"] = ref.get("displayName", "")
            
def run_reactome_hif1_workflow() -> Path:
    """
    Reactome HIF-1 workflow.

    HIF-1 관련 전체 구조를 얻기 위해 다음을 조합한다.

      1. Pathway 중심 실행
         R-HSA-1234174 Cellular response to hypoxia

      2. Entity seed 실행
         R-HSA-1234135 HIF1A [nucleoplasm]

      3. Reaction detail 실행
         R-HSA-1234171 HIF-alpha binds ARNT

      4. Chemical node 추출
         SimpleEntity / ChemicalDrug / ChemicalCompound

    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path("results") / f"reactome_hif1_combined_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    master_graph = nx.DiGraph()
    master_graph.graph["source"] = "reactome"
    master_graph.graph["workflow"] = "reactome_hif1_combined"

    steps = [
        {
            "name": "01_pathway_cellular_response_to_hypoxia",
            "reactome_id": "R-HSA-1234174",
            "description": "Cellular response to hypoxia",
            "build_reactome_flow": True,
            "recursive_events": True,
            "recursive_max_depth": 3,
        },
        {
            "name": "02_entity_hif1a_nucleoplasm",
            "reactome_id": "R-HSA-1234135",
            "description": "HIF1A [nucleoplasm]",
            "build_reactome_flow": True,
            "recursive_events": False,
            "recursive_max_depth": 0,
        },
        {
            "name": "03_reaction_hif_alpha_binds_arnt",
            "reactome_id": "R-HSA-1234171",
            "description": "HIF-alpha binds ARNT forming HIF-alpha:ARNT",
            "build_reactome_flow": True,
            "recursive_events": False,
            "recursive_max_depth": 0,
        },
    ]

    step_summaries = []

    for step_index, step in enumerate(steps, start=1):
        print(f"\n[{step_index}] Fetching {step['reactome_id']} - {step['description']}")

        step_output_dir = output_dir / step["name"]

        try:
            result = fetch_and_merge(
                master_graph,
                step["reactome_id"],
                build_reactome_flow=step["build_reactome_flow"],
                recursive_events=step["recursive_events"],
                recursive_max_depth=step["recursive_max_depth"],
                step_output_dir=step_output_dir,
            )

            step_summaries.append({
                "step": step_index,
                "name": step["name"],
                "reactome_id": step["reactome_id"],
                "description": step["description"],
                "status": result.get("status"),
                "mode": result.get("mode"),
                "node_count": result.get("node_count"),
                "edge_count": result.get("edge_count"),
                "result": result,
            })

        except Exception as exc:
            step_summaries.append({
                "step": step_index,
                "name": step["name"],
                "reactome_id": step["reactome_id"],
                "description": step["description"],
                "status": "failed",
                "error": str(exc),
            })

            print(f"FAILED: {step['reactome_id']} - {exc}")
    
    # 통합 graph 최종 보강
    pathway_for_enrich = SignalingPathway(source="reactome", timeout=(20, 180))
    pathway_for_enrich.graph = master_graph
    enrich_reactome_node_attributes(pathway_for_enrich)
    master_graph = pathway_for_enrich.graph

    # 통합 graph 저장
    save_graph_tables(master_graph, output_dir)

    # relation_type별 결과 저장
    save_relation_edges(master_graph, output_dir)

    # ChemicalCompound / SimpleEntity 추출
    save_node_type_filter(
        master_graph,
        output_dir / "chemical_nodes.csv",
        allowed_types={"SimpleEntity", "ChemicalDrug", "ChemicalCompound"},
    )

    # Complex 추출
    save_node_type_filter(
        master_graph,
        output_dir / "complex_nodes.csv",
        allowed_types={"Complex"},
    )

    # Reaction/Event 추출
    save_node_type_filter(
        master_graph,
        output_dir / "reaction_event_nodes.csv",
        allowed_types={"Reaction", "ReactionLikeEvent", "BlackBoxEvent"},
    )
    
    # pathway 추출
    save_node_type_filter(
        master_graph,
        output_dir / "pathway_nodes.csv",
        allowed_types={"Pathway", "TopLevelPathway"},
    )

    # HIF-1 관련 keyword 검색
    keyword_hits = find_nodes_by_keywords(
        master_graph,
        keywords=[
            "HIF1A",
            "HIF-1",
            "HIF",
            "ARNT",
            "VEGF",
            "VEGFA",
            "EPO",
            "CA9",
            "O2",
            "2OG",
            "CO2",
            "SUCCA",
            "VHL",
            "PHD",
        ],
    )

    save_keyword_hits(
        master_graph,
        keyword_hits,
        output_dir / "reactome_hif1_keyword_hits.csv",
    )

    # 대표 path 확인
    # HIF1A entity -> HIF-alpha binds ARNT reaction
    save_paths(
        master_graph,
        start_node="REACTOME:R-HSA-1234135",
        end_node="REACTOME:R-HSA-1234171",
        output_path=output_dir / "paths_HIF1A_to_HIF_alpha_binds_ARNT_directed.csv",
        cutoff=10,
        undirected=False,
    )

    save_paths(
        master_graph,
        start_node="REACTOME:R-HSA-1234135",
        end_node="REACTOME:R-HSA-1234171",
        output_path=output_dir / "paths_HIF1A_to_HIF_alpha_binds_ARNT_undirected.csv",
        cutoff=10,
        undirected=True,
    )

    summary = {
        "workflow": "reactome_hif1_combined",
        "steps": step_summaries,
        "combined_node_count": master_graph.number_of_nodes(),
        "combined_edge_count": master_graph.number_of_edges(),
        "output_dir": str(output_dir),
        "note": (
            "Reactome HIF-1 workflow combines pathway-centric, entity-seed, "
            "reaction-detail, and chemical-node extraction."
        ),
    }

    with (output_dir / "workflow_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    try:
        pathway_for_plot = SignalingPathway(source="reactome")
        pathway_for_plot.graph = master_graph
        pathway_for_plot.visualize_nodesedges(
            with_edge_labels=False,
            save_path=output_dir / "reactome_hif1_combined_graph.png",
            show=False,
        )
    except Exception as exc:
        print(f"combined graph image skipped: {exc}")

    print(f"\nReactome HIF-1 combined workflow saved: {output_dir}")
    return output_dir


if __name__ == "__main__":
    run_reactome_hif1_workflow()