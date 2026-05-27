# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:30:58 2026

@author: user
"""

# workflows/generic_reactome_target_workflow.py

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


def enrich_reactome_node_attributes(pathway: SignalingPathway) -> None:
    """
    Reactome graph node 중 속성이 부족한 node를 /query/{raw_id}로 다시 조회하여
    schemaClass, compartment, species, reference 정보를 보강한다.
    """
    if pathway.graph.graph.get("source") != "reactome":
        return

    for node_id, attrs in list(pathway.graph.nodes(data=True)):
        raw_id = attrs.get("raw_id")

        if not raw_id:
            continue

        raw_id = str(raw_id)

        if not raw_id.startswith("R-"):
            continue

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

        compartments = detail.get("compartment") or detail.get("compartments") or []
        if isinstance(compartments, dict):
            compartments = [compartments]

        for comp in compartments:
            if isinstance(comp, dict):
                attrs["compartment"] = comp.get("displayName", comp.get("name", ""))
                break

        species = detail.get("species")
        if isinstance(species, dict):
            attrs["species"] = species.get("displayName", species.get("name", ""))

        ref = detail.get("referenceEntity")
        if isinstance(ref, dict):
            attrs["reference_identifier"] = ref.get("identifier", "")
            attrs["reference_database"] = ref.get("databaseName", "")
            attrs["reference_display"] = ref.get("displayName", "")


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
        "entity_complex_links.csv": {
            "component_of",
            "component",
            "member",
            "candidate",
        },
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


def find_nodes_by_keywords(
    graph: nx.DiGraph,
    keywords: list[str],
) -> dict[str, list[str]]:
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


def save_keyword_hits(
    graph: nx.DiGraph,
    keyword_hits: dict[str, list[str]],
    output_path: Path,
) -> None:
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
        paths = list(
            nx.all_simple_paths(
                search_graph,
                start_node,
                end_node,
                cutoff=cutoff,
            )
        )
    else:
        paths = []

    with output_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "path_index",
            "step_index",
            "node_id",
            "label",
            "type",
            "relation_to_next",
        ])

        for path_index, path in enumerate(paths, start=1):
            for step_index, node_id in enumerate(path, start=1):
                relation_to_next = ""

                if step_index < len(path):
                    next_node = path[step_index]

                    if graph.has_edge(node_id, next_node):
                        relation_to_next = graph.edges[node_id, next_node].get(
                            "relation_type",
                            "",
                        )

                    elif graph.has_edge(next_node, node_id):
                        relation_to_next = (
                            "reverse:"
                            + graph.edges[next_node, node_id].get("relation_type", "")
                        )

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
    step: dict,
    *,
    step_output_dir: Path,
) -> dict:
    """
    하나의 Reactome step을 실행하고 master_graph에 병합한다.
    """
    reactome_id = step["reactome_id"]

    pathway = SignalingPathway(source="reactome", timeout=(20, 180))

    result = pathway.fetch_pathway(
        reactome_id,
        build_reactome_flow=step.get("build_reactome_flow", True),
        recursive_events=step.get("recursive_events", False),
        recursive_max_depth=step.get("recursive_max_depth", 3),
    )

    # 핵심: 저장 전에 node 속성 보강
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

    master_graph.update(pathway.graph)

    return result


def run_generic_reactome_target_workflow(
    target_key: str,
    target_config: dict,
) -> Path | None:
    """
    target_registry.py의 reactome.steps를 이용하여
    target별 Reactome workflow를 일반 실행한다.
    """
    reactome_config = target_config.get("reactome", {})
    steps = reactome_config.get("steps", [])

    if not steps:
        print(f"\n{target_key}에 대한 Reactome workflow steps가 등록되어 있지 않습니다.")
        print("target_registry.py의 reactome.steps를 먼저 추가해야 합니다.")
        return None

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path("results") / f"reactome_{target_key}_combined_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    master_graph = nx.DiGraph()
    master_graph.graph["source"] = "reactome"
    master_graph.graph["workflow"] = f"reactome_{target_key}_combined"

    step_summaries = []

    for step_index, step in enumerate(steps, start=1):
        print(
            f"\n[{step_index}] Fetching {step['reactome_id']} - "
            f"{step.get('description', '')}"
        )

        step_name = step.get("name", f"step_{step_index:02d}_{step['reactome_id']}")
        safe_step_name = (
            step_name
            .replace(":", "_")
            .replace("/", "_")
            .replace("\\", "_")
        )

        step_output_dir = output_dir / safe_step_name

        try:
            result = fetch_and_merge(
                master_graph,
                step,
                step_output_dir=step_output_dir,
            )

            step_summaries.append({
                "step": step_index,
                "name": step_name,
                "reactome_id": step["reactome_id"],
                "description": step.get("description", ""),
                "status": result.get("status"),
                "mode": result.get("mode"),
                "node_count": result.get("node_count"),
                "edge_count": result.get("edge_count"),
                "result": result,
            })

        except Exception as exc:
            step_summaries.append({
                "step": step_index,
                "name": step_name,
                "reactome_id": step.get("reactome_id", ""),
                "description": step.get("description", ""),
                "status": "failed",
                "error": str(exc),
            })

            print(f"FAILED: {step.get('reactome_id', '')} - {exc}")

    # 통합 graph 최종 보강
    pathway_for_enrich = SignalingPathway(source="reactome", timeout=(20, 180))
    pathway_for_enrich.graph = master_graph
    enrich_reactome_node_attributes(pathway_for_enrich)
    master_graph = pathway_for_enrich.graph

    # 통합 graph 저장
    save_graph_tables(master_graph, output_dir)

    # relation_type별 분리 저장
    save_relation_edges(master_graph, output_dir)

    # node type별 분리 저장
    save_node_type_filter(
        master_graph,
        output_dir / "chemical_nodes.csv",
        allowed_types={"SimpleEntity", "ChemicalDrug", "ChemicalCompound"},
    )

    save_node_type_filter(
        master_graph,
        output_dir / "complex_nodes.csv",
        allowed_types={"Complex"},
    )

    save_node_type_filter(
        master_graph,
        output_dir / "reaction_event_nodes.csv",
        allowed_types={"Reaction", "ReactionLikeEvent", "BlackBoxEvent"},
    )

    save_node_type_filter(
        master_graph,
        output_dir / "pathway_nodes.csv",
        allowed_types={"Pathway", "TopLevelPathway"},
    )

    # keyword hits
    keyword_hits = find_nodes_by_keywords(
        master_graph,
        reactome_config.get("keywords", []),
    )

    save_keyword_hits(
        master_graph,
        keyword_hits,
        output_dir / f"reactome_{target_key}_keyword_hits.csv",
    )

    # representative paths
    path_summaries = []

    for path_cfg in reactome_config.get("representative_paths", []):
        name = path_cfg["name"]
        start = path_cfg["start"]
        end = path_cfg["end"]

        directed_paths = save_paths(
            master_graph,
            start_node=start,
            end_node=end,
            output_path=output_dir / f"paths_directed_{name}.csv",
            cutoff=10,
            undirected=False,
        )

        undirected_paths = save_paths(
            master_graph,
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
        "workflow": "generic_reactome_target",
        "target_key": target_key,
        "steps": step_summaries,
        "path_summaries": path_summaries,
        "combined_node_count": master_graph.number_of_nodes(),
        "combined_edge_count": master_graph.number_of_edges(),
        "output_dir": str(output_dir),
        "note": (
            "Generic Reactome workflow combines pathway-centric, entity-seed, "
            "reaction-detail, and chemical-node extraction based on target_registry.py."
        ),
    }

    with (output_dir / "workflow_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    try:
        pathway_for_plot = SignalingPathway(source="reactome")
        pathway_for_plot.graph = master_graph
        pathway_for_plot.visualize_nodesedges(
            with_edge_labels=False,
            save_path=output_dir / f"reactome_{target_key}_combined_graph.png",
            show=False,
        )
    except Exception as exc:
        print(f"combined graph image skipped: {exc}")

    print(f"\nReactome {target_key} workflow saved: {output_dir}")
    return output_dir