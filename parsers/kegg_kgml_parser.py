"""KEGG KGML parser.

Converts KEGG KGML XML into a directed NetworkX graph. Node IDs are prefixed
with ``KEGG:`` so they cannot be confused with Reactome stable IDs.
"""

from __future__ import annotations

from typing import Any
from lxml import etree
import networkx as nx
#from bs4 import BeautifulSoup

def _expand_group_nodes(
    node_id: str,
    graph: nx.DiGraph,
) -> list[str]:
    """
    group node이면 group_components를 풀어서 component node list를 반환.
    group node가 아니면 자기 자신만 반환.
    """
    if node_id not in graph:
        return [node_id]

    attrs = graph.nodes[node_id]
    if attrs.get("type") != "group":
        return [node_id]

    components = str(attrs.get("group_components", "")).split(";")
    components = [x for x in components if x]

    return components or [node_id]

def _clean_kegg_name(name: str | None) -> list[str]:
    if not name:
        return []

    # 예: "hsa:3091 hsa:405"
    tokens = str(name).strip().split()
    return tokens


def _kegg_node_id_from_name(name: str, entry_type: str, pathway_id: str, entry_id: str) -> str:
    """
    KGML entry name을 graph node_id로 변환한다.
    """
    if entry_type == "gene":
        # hsa:3091 -> KEGG:3091
        first = _clean_kegg_name(name)[0] if _clean_kegg_name(name) else entry_id
        return f"KEGG:{first.split(':')[-1]}"

    if entry_type == "compound":
        # cpd:C00007 또는 C00007 -> KEGG:C00007
        first = _clean_kegg_name(name)[0] if _clean_kegg_name(name) else entry_id
        return f"KEGG:{first.split(':')[-1]}"

    if entry_type == "map":
        first = _clean_kegg_name(name)[0] if _clean_kegg_name(name) else entry_id
        return f"KEGG:{first.split(':')[-1]}"

    if entry_type == "group":
        return f"KEGG_GROUP:{pathway_id}:{entry_id}"

    return f"KEGG_ENTRY:{pathway_id}:{entry_id}"

def _get_entry_label(entry: etree._Element) -> str:
    graphics = entry.find("graphics")

    if graphics is not None:
        label = graphics.get("name")
        if label:
            return label

    name = entry.get("name")
    if name:
        return name

    return entry.get("id", "unknown")


def _build_group_label(
    graph: nx.DiGraph,
    component_node_ids: list[str],
    fallback: str,
) -> str:
    labels = []

    for node_id in component_node_ids:
        if node_id in graph:
            label = graph.nodes[node_id].get("label", node_id)
            labels.append(str(label))

    if labels:
        return "Group: " + " / ".join(labels[:5])

    return fallback if fallback and fallback != "undefined" else "Group"


def parse_kegg_kgml(
    kgml_data: str | bytes,
    graph: nx.DiGraph,
    pathway_id: str,
    expand_group_relations: bool = True,
) -> dict[str, Any]:
    """
    KEGG KGML을 graph로 변환한다.

    - group entry를 KEGG_GROUP:{pathway_id}:{entry_id} node로 보존
    - group_components / group_component_count를 node attribute로 저장
    - group -> component edge를 group_component relation으로 저장
    - expand_group_relations=True이면 GROUP이 관여한 relation을 실제 component node까지 via_group edge로 확장
    """
    added_edges = 0
    group_component_edge_count = 0
    group_expanded_edge_count = 0

    if isinstance(kgml_data, str):
        root = etree.fromstring(kgml_data.encode("utf-8"))
    else:
        root = etree.fromstring(kgml_data)

    graph.graph.update({
        "source": "kegg",
        "pathway_id": pathway_id,
        "name": root.get("title", pathway_id),
    })

    entries = root.findall("entry")
    relations = root.findall("relation")

    entry_id_to_node_id: dict[str, str] = {}
    entry_id_to_entry: dict[str, etree._Element] = {}
    group_components: dict[str, list[str]] = {}

    # ------------------------------------------------------------
    # 1st pass: group이 아닌 node 먼저 추가
    # ------------------------------------------------------------
    for entry in entries:
        entry_id = entry.get("id")
        entry_type = entry.get("type", "")
        entry_name = entry.get("name", "")

        if not entry_id:
            continue

        entry_id_to_entry[entry_id] = entry

        if entry_type == "group":
            component_ids = [
                comp.get("id")
                for comp in entry.findall("component")
                if comp.get("id")
            ]
            group_components[entry_id] = component_ids
            continue

        node_id = _kegg_node_id_from_name(
            entry_name,
            entry_type,
            pathway_id,
            entry_id,
        )

        label = _get_entry_label(entry)

        graph.add_node(
            node_id,
            source="kegg",
            type=entry_type,
            label=label,
            raw_id=node_id.replace("KEGG:", ""),
            kgml_entry_id=entry_id,
            kgml_name=entry_name,
        )

        entry_id_to_node_id[entry_id] = node_id

    # ------------------------------------------------------------
    # 2nd pass: group node 추가 + component edge 추가
    # ------------------------------------------------------------
    for entry_id, component_entry_ids in group_components.items():
        entry = entry_id_to_entry[entry_id]

        group_node_id = _kegg_node_id_from_name(
            entry.get("name", ""),
            "group",
            pathway_id,
            entry_id,
        )

        component_node_ids = [
            entry_id_to_node_id[component_entry_id]
            for component_entry_id in component_entry_ids
            if component_entry_id in entry_id_to_node_id
        ]

        fallback_label = _get_entry_label(entry)
        group_label = _build_group_label(
            graph,
            component_node_ids,
            fallback=fallback_label,
        )

        graph.add_node(
            group_node_id,
            source="kegg",
            type="group",
            label=group_label,
            raw_id=entry_id,
            kgml_entry_id=entry_id,
            kgml_name=entry.get("name", ""),
            group_components=";".join(component_node_ids),
            group_component_count=len(component_node_ids),
        )

        entry_id_to_node_id[entry_id] = group_node_id

        for component_node_id in component_node_ids:
            graph.add_edge(
                group_node_id,
                component_node_id,
                relation_type="group_component",
                source="kegg",
            )
            group_component_edge_count += 1

    # ------------------------------------------------------------
    # 3rd pass: relation 추가
    # ------------------------------------------------------------
    for relation in relations:
        entry1 = relation.get("entry1")
        entry2 = relation.get("entry2")

        if not entry1 or not entry2:
            continue

        source_node = entry_id_to_node_id.get(entry1)
        target_node = entry_id_to_node_id.get(entry2)

        if not source_node or not target_node:
            continue

        subtype_values = []
        for subtype in relation.findall("subtype"):
            subtype_name = subtype.get("name")
            if subtype_name:
                subtype_values.append(subtype_name)

        relation_type = ", ".join(subtype_values) if subtype_values else relation.get("type", "")

        graph.add_edge(
            source_node,
            target_node,
            relation_type=relation_type,
            source="kegg",
            kgml_relation_type=relation.get("type", ""),
        )
        added_edges += 1

        if not expand_group_relations:
            continue

        source_expanded = _expand_group_nodes(source_node, graph)
        target_expanded = _expand_group_nodes(target_node, graph)

        for expanded_source in source_expanded:
            for expanded_target in target_expanded:
                if expanded_source == source_node and expanded_target == target_node:
                    continue

                if expanded_source == expanded_target:
                    continue

                graph.add_edge(
                    expanded_source,
                    expanded_target,
                    relation_type=f"{relation_type}, via_group",
                    source="kegg",
                    kgml_relation_type=relation.get("type", ""),
                    expanded_from_group=True,
                    original_source_node=source_node,
                    original_target_node=target_node,
                )
                group_expanded_edge_count += 1

    return {
        "status": "success",
        "source": "kegg",
        "pathway_id": pathway_id,
        "node_count": len(graph.nodes),
        "edge_count": len(graph.edges),
        "added_edges": added_edges,
        "group_count": len(group_components),
        "group_component_edge_count": group_component_edge_count,
        "group_expanded_edge_count": group_expanded_edge_count,
        "expand_group_relations": expand_group_relations,
    }
