"""Reactome JSON parser.

Converts Reactome ContentService JSON into a directed NetworkX graph. Node IDs
are prefixed with ``REACTOME:`` so they cannot be confused with KEGG IDs.

Created on Mon Sep  8 14:15:50 2025

@author: USER, SANG JIN PARK
"""

from __future__ import annotations

from typing import Any

import networkx as nx


def _iter_dicts(value: Any) -> list[dict[str, Any]]:
    """
    Reactome 응답에서 dict 항목만 안전하게 순회하기 위한 helper.

    허용:
      None        -> []
      dict        -> [dict]
      list[dict]  -> dict만 필터링
      int/str/etc -> []

    Reactome API 응답에는 경우에 따라 int dbId, 문자열, None 등이 섞일 수 있으므로
    parser에서는 반드시 이 함수를 통해 순회하는 것이 안전하다.
    """
    if value is None:
        return []

    if isinstance(value, dict):
        return [value]

    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]

    return []


def _reactome_node_id(stable_id: str | None) -> str | None:
    if not stable_id:
        return None

    stable_id = str(stable_id)

    if stable_id.startswith("REACTOME:"):
        return stable_id

    return f"REACTOME:{stable_id}"

def _first_compartment_name(entity: dict[str, Any]) -> str:
    compartments = entity.get("compartment") or entity.get("compartments")
    for comp in _iter_dicts(compartments):
        name = comp.get("displayName") or comp.get("name")
        if name:
            return str(name)
    return ""

def _species_name(entity: dict[str, Any]) -> str:
    for species in _iter_dicts(entity.get("species")):
        name = species.get("displayName") or species.get("name")
        if name:
            return str(name)
    return ""

def _reference_info(entity: dict[str, Any]) -> tuple[str, str, str]:
    refs = _iter_dicts(entity.get("referenceEntity")) or _iter_dicts(entity.get("referenceEntities"))
    if not refs:
        return "", "", ""
    ref = refs[0]
    identifier = str(ref.get("identifier", ""))
    database = str(ref.get("databaseName", ""))
    display = str(ref.get("displayName", ""))
    return identifier, database, display

def _add_reactome_entity_node(
    graph: nx.DiGraph,
    entity: dict[str, Any] | None,
) -> str | None:
    """
    Reactome entity/event/reaction node를 graph에 추가.

    entity가 dict가 아니면 None 반환.
    이 방어가 없으면 int object has no attribute get 에러가 발생할 수 있다.
    """
    if not isinstance(entity, dict):
        return None

    raw_id = entity.get("stId") or entity.get("dbId")
    if raw_id is None:
        return None

    raw_id = str(raw_id)
    node_id = _reactome_node_id(raw_id)

    if node_id is None:
        return None

    ref_id, ref_db, ref_display = _reference_info(entity)

    graph.add_node(
        node_id,
        label=entity.get("displayName", raw_id),
        type=entity.get("schemaClass", "Unknown"),
        source="reactome",
        raw_id=raw_id,
        schemaClass=entity.get("schemaClass", "Unknown"),
        compartment=_first_compartment_name(entity),
        species=_species_name(entity),
        reference_identifier=ref_id,
        reference_database=ref_db,
        reference_display=ref_display,
    )

    return node_id


def parse_reactome_events(
    events: list[dict[str, Any]] | dict[str, Any] | Any,
    graph: nx.DiGraph,
    pathway_id: str | None = None,
) -> dict[str, Any]:
    """Parse containedEvents or a single-query payload into graph nodes.

    This stage builds event/entity nodes. Detailed input/output/catalyst edges
    are added by ``parse_reactome_reaction_detail`` when detailed reaction data
    are available.
    """
    graph.graph.update(
        {
            "source": "reactome",
            "pathway_id": pathway_id,
            "name": pathway_id or "Reactome Pathway",
        }
    )

    parsed_nodes: list[dict[str, Any]] = []
    edge_count_before = len(graph.edges)

    for entity in _iter_dicts(events):
        node_id = _add_reactome_entity_node(graph, entity)
        if not node_id:
            continue

        parsed_nodes.append(
            {
                "id": node_id,
                "raw_id": graph.nodes[node_id].get("raw_id"),
                "label": graph.nodes[node_id].get("label"),
                "type": graph.nodes[node_id].get("type"),
                "source": "reactome",
            }
        )

        # Some Reactome objects contain refEntities or hasComponent-like fields.
        for field_name, rel_type in (
            ("refEntities", "reference"),
            ("hasComponent", "component"),
            ("hasMember", "member"),
            ("hasCandidate", "candidate"),
        ):
            for child in _iter_dicts(entity.get(field_name)):
                child_id = _add_reactome_entity_node(graph, child)
                if child_id:
                    graph.add_edge(
                        node_id,
                        child_id,
                        relation_type=rel_type,
                        source="reactome",
                    )

    return {
        "status": "success",
        "source": "reactome",
        "pathway_id": pathway_id,
        "nodes": parsed_nodes,
        "added_nodes": len(parsed_nodes),
        "added_edges": len(graph.edges) - edge_count_before,
        "node_count": len(graph.nodes),
        "edge_count": len(graph.edges),
    }


def parse_reactome_reaction_detail(
    reaction_detail: dict[str, Any] | Any,
    graph: nx.DiGraph,
) -> dict[str, Any]:
    """Parse a detailed Reactome ReactionLikeEvent into input/output/catalyst edges."""

    if not isinstance(reaction_detail, dict):
        return {
            "status": "skipped",
            "reason": "reaction_detail is not a dict",
            "source": "reactome",
            "added_edges": 0,
            "node_count": len(graph.nodes),
            "edge_count": len(graph.edges),
        }

    reaction_id = _add_reactome_entity_node(graph, reaction_detail)
    if not reaction_id:
        return {
            "status": "skipped",
            "reason": "reaction has no stId/dbId",
            "source": "reactome",
            "added_edges": 0,
            "node_count": len(graph.nodes),
            "edge_count": len(graph.edges),
        }

    graph.nodes[reaction_id]["type"] = reaction_detail.get("schemaClass", "ReactionLikeEvent")
    edge_count_before = len(graph.edges)

    for input_entity in _iter_dicts(reaction_detail.get("input")):
        input_id = _add_reactome_entity_node(graph, input_entity)
        if input_id:
            graph.add_edge(input_id, reaction_id, relation_type="input", source="reactome")

    for output_entity in _iter_dicts(reaction_detail.get("output")):
        output_id = _add_reactome_entity_node(graph, output_entity)
        if output_id:
            graph.add_edge(reaction_id, output_id, relation_type="output", source="reactome")

    for catalyst in _iter_dicts(reaction_detail.get("catalystActivity")):
        for physical_entity in _iter_dicts(catalyst.get("physicalEntity")):
            catalyst_id = _add_reactome_entity_node(graph, physical_entity)
            if catalyst_id:
                graph.add_edge(catalyst_id, reaction_id, relation_type="catalyst", source="reactome")

    for regulation in _iter_dicts(reaction_detail.get("regulatedBy")):
        relation_type = regulation.get("schemaClass", "regulator")
        for regulator in _iter_dicts(regulation.get("regulator")):
            regulator_id = _add_reactome_entity_node(graph, regulator)
            if regulator_id:
                graph.add_edge(regulator_id, reaction_id, relation_type=relation_type, source="reactome")

    # Event ordering: precedingEvent -> reaction_id
    for preceding in _iter_dicts(reaction_detail.get("precedingEvent")):
        preceding_id = _add_reactome_entity_node(graph, preceding)
        if preceding_id:
            graph.add_edge(preceding_id, reaction_id, relation_type="precedingEvent", source="reactome")

    return {
        "status": "success",
        "source": "reactome",
        "reaction_id": reaction_id,
        "added_edges": len(graph.edges) - edge_count_before,
        "node_count": len(graph.nodes),
        "edge_count": len(graph.edges),
    }
