"""
Core pathway manager.

The public class ``SignalingPathway`` provides one source-aware entry point:
``fetch_pathway``. KEGG and Reactome-specific network calls and parsers are
located in separate provider/parser modules.

Created on Mon Sep  8 14:15:50 2025

@author: USER, SANG JIN PARK
"""

from __future__ import annotations

import json, csv
from enum import Enum
from pathlib import Path
from typing import Any
from datetime import datetime

import networkx as nx

from parsers.kegg_kgml_parser import parse_kegg_kgml
from parsers.reactome_json_parser import (
    parse_reactome_events,
    parse_reactome_reaction_detail,
    _iter_dicts,
)

from providers.kegg_provider import KeggProvider
from providers.reactome_provider import ReactomeProvider
from visualization.pathway_visualizer import visualize_graph


class PathwaySource(str, Enum):
    KEGG = "kegg"
    REACTOME = "reactome"


class SignalingPathway:
    """Source-aware pathway graph manager."""

    def __init__(self, source: str = "kegg", timeout: int = 30):
        self.source = self._normalize_source(source)
        self.graph = nx.DiGraph()
        self.graph.graph["source"] = self.source.value
        self.kegg_provider = KeggProvider(timeout=timeout)
        self.reactome_provider = ReactomeProvider(timeout=timeout)

    @staticmethod
    def _normalize_source(source: str) -> PathwaySource:
        source_norm = str(source).strip().lower()
        if source_norm in ("kegg", "kgml"):
            return PathwaySource.KEGG
        if source_norm in ("reactome", "react"):
            return PathwaySource.REACTOME
        raise ValueError("source는 'kegg' 또는 'reactome'만 허용됩니다.")

    def normalize_node_id(self, node_id):
        """
        사용자가 3091처럼 입력해도 KEGG:3091로 자동 변환.
        Reactome ID도 REACTOME: prefix를 자동 보정.
        """
        node_id = str(node_id).strip()
    
        if node_id in self.graph:
            return node_id
    
        source = self.graph.graph.get(
            "source", 
            self.source.value if hasattr(self.source, "value") else self.source
        )
        source = str(source).lower()
    
        candidates = []
    
        if source == "kegg":
            candidates = [
                f"KEGG:{node_id}",
                node_id.replace("hsa:", "KEGG:"),
                node_id.replace("path:", "KEGG:")
            ]
    
        elif source == "reactome":
            candidates = [
                f"REACTOME:{node_id}",
                node_id.replace("reactome:", "REACTOME:")
            ]
    
        for candidate in candidates:
            if candidate in self.graph:
                return candidate
    
        return node_id

    def validate_pathway_id(self, pathway_id: str) -> str:
        """Validate pathway ID based on current data source."""
        if self.source == PathwaySource.KEGG:
            return self.kegg_provider.normalize_pathway_id(pathway_id)
        if self.source == PathwaySource.REACTOME:
            return self.reactome_provider.normalize_pathway_id(pathway_id)
        raise ValueError(f"Unsupported source: {self.source}")

    def reset_graph(self, pathway_id: str | None = None) -> None:
        self.graph.clear()
        self.graph.graph["source"] = self.source.value
        if pathway_id:
            self.graph.graph["pathway_id"] = pathway_id

    def fetch_reactome_pathway_flow(
        self,
        pathway_id: str,
        pathway_payload: dict[str, Any] | None = None,
        build_reactome_flow: bool = True,
        recursive_events: bool = False,
        recursive_max_depth: int = 3,
    ) -> dict[str, Any]:
        """
        Reactome Pathway ID를 받아 containedEvents를 graph로 구성한다.
    
        recursive_events=False:
          - 현재 pathway의 immediate containedEvents만 사용
    
        recursive_events=True:
          - 하위 Pathway의 containedEvents까지 재귀적으로 확장
          - HIF-1처럼 하위 pathway/event가 깊은 경우 유용
        """
        if self.source != PathwaySource.REACTOME:
            raise ValueError(
                "fetch_reactome_pathway_flow()는 source='reactome'일 때만 호출할 수 있습니다."
            )
    
        clean_id = self.validate_pathway_id(pathway_id)
    
        if pathway_payload is None:
            pathway_payload = self.reactome_provider.query_reactome_id(clean_id)
    
        payload = self.reactome_provider.fetch_pathway_payload(clean_id)
    
        if recursive_events:
            raw_events = self.reactome_provider.fetch_contained_events_recursive(
                clean_id,
                max_depth=recursive_max_depth,
            )
        else:
            raw_events = payload.get("events", [])
    
        events = _iter_dicts(raw_events)
    
        pathway_node_id = f"REACTOME:{clean_id}"
    
        self.graph.add_node(
            pathway_node_id,
            label=pathway_payload.get("displayName", clean_id),
            type=pathway_payload.get("schemaClass", "Pathway"),
            source="reactome",
            raw_id=clean_id,
            schemaClass=pathway_payload.get("schemaClass", "Pathway"),
        )
    
        added_edges = 0
        added_reactions = 0
        skipped_non_dict_events = len(raw_events) - len(events)
    
        summary = parse_reactome_events(events, self.graph, clean_id)
    
        for event in events:
            event_id = event.get("stId")
            if not event_id:
                continue
    
            schema_class = event.get("schemaClass", "")
            event_node_id = f"REACTOME:{event_id}"
    
            if event_node_id not in self.graph:
                self.graph.add_node(
                    event_node_id,
                    label=event.get("displayName", event_id),
                    type=schema_class,
                    source="reactome",
                    raw_id=event_id,
                    schemaClass=schema_class,
                )
    
            self.graph.add_edge(
                pathway_node_id,
                event_node_id,
                relation_type="has_event",
                source="reactome",
            )
            added_edges += 1
    
            if not build_reactome_flow:
                continue
    
            if schema_class not in {
                "ReactionLikeEvent",
                "Reaction",
                "BlackBoxEvent",
                "Polymerisation",
                "Depolymerisation",
            }:
                continue
    
            detail = self.reactome_provider.fetch_reaction_detail(event_id)
    
            if detail:
                result = parse_reactome_reaction_detail(detail, self.graph)
                added_edges += result.get("added_edges", 0)
                added_reactions += 1
    
        summary.update(
            {
                "status": "success",
                "source": "reactome",
                "mode": "pathway_flow_recursive" if recursive_events else "pathway_flow",
                "pathway_id": clean_id,
                "schemaClass": pathway_payload.get("schemaClass", "Pathway"),
                "displayName": pathway_payload.get("displayName", clean_id),
                "recursive_events": recursive_events,
                "recursive_max_depth": recursive_max_depth if recursive_events else 0,
                "raw_event_count": len(raw_events),
                "event_count": len(events),
                "skipped_non_dict_events": skipped_non_dict_events,
                "added_reactions": added_reactions,
                "added_edges": added_edges,
                "node_count": len(self.graph.nodes),
                "edge_count": len(self.graph.edges),
            }
        )
    
        return summary

    def fetch_reactome(
        self,
        reactome_id: str,
        build_reactome_flow: bool = True,
        recursive_events: bool = False,
        recursive_max_depth: int = 3,
    ) -> dict[str, Any]:
        """
        Reactome ID의 schemaClass를 확인한 뒤,
        Pathway / Reaction / Entity 처리 함수로 분기한다.
    
        TopLevelPathway는 범위가 매우 크므로 기본적으로 reaction detail 확장을 끄고,
        Pathway → Event 계층 중심으로만 구성한다.
        """
        if self.source != PathwaySource.REACTOME:
            raise ValueError("fetch_reactome()는 source='reactome'일 때만 호출할 수 있습니다.")
    
        clean_id = self.validate_pathway_id(reactome_id)
        data = self.reactome_provider.query_reactome_id(clean_id)
        schema_class = data.get("schemaClass")
    
        if not schema_class:
            raise ValueError(f"Reactome schemaClass를 확인할 수 없습니다: {clean_id}")
    
        self.graph.graph["source"] = "reactome"
        self.graph.graph["pathway_id"] = clean_id
        self.graph.graph["reactome_schemaClass"] = schema_class
        self.graph.graph["name"] = data.get("displayName", clean_id)
    
        pathway_schema_classes = {
            "Pathway",
            "TopLevelPathway",
        }
    
        entity_schema_classes = {
            "EntityWithAccessionedSequence",
            "Complex",
            "CandidateSet",
            "DefinedSet",
            "OpenSet",
            "GenomeEncodedEntity",
            "ChemicalDrug",
            "SimpleEntity",
            "Polymer",
        }
    
        reaction_schema_classes = {
            "Reaction",
            "ReactionLikeEvent",
            "BlackBoxEvent",
            "Polymerisation",
            "Depolymerisation",
            "FailedReaction",
        }
    
        if schema_class in pathway_schema_classes:
            if schema_class == "TopLevelPathway":
                print(
                    "ℹ️ TopLevelPathway입니다. "
                    "범위가 크므로 reaction detail 확장은 생략하고 "
                    "Pathway → Event 계층만 우선 구성합니다."
                )
                build_reactome_flow = False
                recursive_events = False
        
            return self.fetch_reactome_pathway_flow(
                clean_id,
                pathway_payload=data,
                build_reactome_flow=build_reactome_flow,
                recursive_events=recursive_events,
                recursive_max_depth=recursive_max_depth,
            )
    
        if schema_class in entity_schema_classes:
            return self.fetch_reactome_entity_neighborhood(
                clean_id,
                seed_payload=data,
                build_reaction_details=build_reactome_flow,
            )
    
        if schema_class in reaction_schema_classes:
            return self.fetch_reactome_reaction_detail(
                clean_id,
                reaction_payload=data,
            )
    
        raise ValueError(f"Unsupported Reactome schemaClass: {schema_class}")

    def enrich_kegg_compound_names(self) -> dict[str, str]:
        """
        KEGG compound node의 label을 KEGG compound name으로 보강한다.
    
        기존:
          KEGG:C00007 label=C00007
    
        보강 후:
          KEGG:C00007 label=Oxygen
          compound_id=C00007
          compound_name=Oxygen
        """
        if self.source != PathwaySource.KEGG:
            return {}
    
        compound_name_map: dict[str, str] = {}
    
        for node_id, attrs in list(self.graph.nodes(data=True)):
            if attrs.get("type") != "compound":
                continue
    
            compound_id_candidate = (
                attrs.get("compound_id")
                or attrs.get("raw_id")
                or str(node_id).replace("KEGG:", "")
            )
    
            try:
                compound_id = self.kegg_provider.normalize_compound_id(
                    str(compound_id_candidate)
                )
            except Exception:
                continue
    
            compound_name = self.kegg_provider.fetch_compound_name(compound_id)
    
            attrs["compound_id"] = compound_id
    
            if not compound_name:
                continue
    
            attrs["compound_name"] = compound_name
            attrs["label"] = compound_name
    
            compound_name_map[compound_id] = compound_name
    
        return compound_name_map

    def fetch_pathway(
        self,
        pathway_id: str,
        build_reactome_flow: bool = True,
        recursive_events: bool = False,
        recursive_max_depth: int = 3,
    ) -> dict[str, Any]:
        """Fetch and parse a pathway/entity/reaction using the configured source.
    
        KEGG:
          - KGML pathway만 처리
    
        Reactome:
          - schemaClass를 먼저 확인한 뒤
            Pathway / ReactionLikeEvent / PhysicalEntity로 자동 분기
        """
        clean_id = self.validate_pathway_id(pathway_id)
        self.reset_graph(clean_id)
    
        if self.source == PathwaySource.KEGG:
            payload = self.kegg_provider.fetch_pathway_payload(clean_id)
            result = parse_kegg_kgml(payload["data"], self.graph, clean_id)
        
            self.enrich_kegg_compound_names()
        
            result.update({
                "node_count": len(self.graph.nodes),
                "edge_count": len(self.graph.edges),
                "compound_name_enriched": True,
            })
        
            return result
    
        if self.source == PathwaySource.REACTOME:
            return self.fetch_reactome(
                clean_id,
                build_reactome_flow=build_reactome_flow,
                recursive_events=recursive_events,
                recursive_max_depth=recursive_max_depth,
            )
    
        raise ValueError(f"Unsupported source: {self.source}")


    # ------------------------------------------------------------------
    # Compatibility wrappers for the older script interface
    # ------------------------------------------------------------------
    def search_pathway_id(self, keyword: str) -> list[dict[str, str]]:
        """Search pathway IDs.

        Currently implemented for KEGG only, matching the original script.
        Reactome search can be added later as a separate provider method.
        """
        if self.source != PathwaySource.KEGG:
            raise NotImplementedError("Reactome keyword search is not implemented in this refactor.")
        return self.kegg_provider.search_pathway_id(keyword)

    def fetch_kegg_pathway(self, pathway_id: str) -> dict[str, Any]:
        """Compatibility wrapper; only valid when source is KEGG."""
        if self.source != PathwaySource.KEGG:
            raise ValueError("fetch_kegg_pathway()는 source='kegg'일 때만 호출할 수 있습니다.")
        return self.fetch_pathway(pathway_id)
    

#    def fetch_reactome_flow(self, pathway_id: str) -> dict[str, Any]:
#        """Compatibility wrapper for Reactome pathway flow."""
#        if self.source != PathwaySource.REACTOME:
#            raise ValueError("fetch_reactome_flow()는 source='reactome'일 때만 호출할 수 있습니다.")
#    
#        return self.fetch_reactome_pathway_flow(
#            pathway_id,
#            build_reactome_flow=True,
#        )


    def fetch_reactome_flow(self, pathway_id: str) -> dict[str, Any]:
        """Compatibility wrapper for Reactome pathway flow."""
        if self.source != PathwaySource.REACTOME:
            raise ValueError("fetch_reactome_flow()는 source='reactome'일 때만 호출할 수 있습니다.")
    
        return self.fetch_reactome_pathway_flow(
            pathway_id,
            build_reactome_flow=True,
        )
    

    def fetch_reactome_pathway(self, pathway_id: str) -> dict[str, Any]:
        """Compatibility wrapper; only valid when source is Reactome."""
        if self.source != PathwaySource.REACTOME:
            raise ValueError("fetch_reactome_pathway()는 source='reactome'일 때만 호출할 수 있습니다.")
    
        return self.fetch_reactome_pathway_flow(
            pathway_id,
            build_reactome_flow=False,
        )
    
    def fetch_reactome_reaction_detail(
        self,
        reaction_id: str,
        reaction_payload: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """
        단일 Reactome ReactionLikeEvent를 중심으로
        input/output/catalyst/regulator graph를 만든다.
        """
        if self.source != PathwaySource.REACTOME:
            raise ValueError("fetch_reactome_reaction_detail()는 source='reactome'일 때만 호출할 수 있습니다.")
    
        clean_id = self.validate_pathway_id(reaction_id)
    
        if reaction_payload is None:
            reaction_payload = self.reactome_provider.fetch_reaction_detail(clean_id)
    
        if not reaction_payload:
            raise ValueError(f"Reactome reaction detail을 가져오지 못했습니다: {clean_id}")
    
        result = parse_reactome_reaction_detail(reaction_payload, self.graph)
    
        reaction_node_id = f"REACTOME:{clean_id}"
        if reaction_node_id in self.graph:
            self.graph.nodes[reaction_node_id]["source"] = "reactome"
            self.graph.nodes[reaction_node_id]["raw_id"] = clean_id
    
        result.update(
            {
                "status": "success",
                "source": "reactome",
                "mode": "reaction_detail",
                "reaction_id": clean_id,
                "schemaClass": reaction_payload.get("schemaClass"),
                "displayName": reaction_payload.get("displayName", clean_id),
                "node_count": len(self.graph.nodes),
                "edge_count": len(self.graph.edges),
            }
        )
    
        return result
    
    def fetch_reactome_entity_neighborhood(
        self,
        entity_id: str,
        seed_payload: dict[str, Any] | None = None,
        build_reaction_details: bool = True,
    ) -> dict[str, Any]:
        """
        Reactome PhysicalEntity를 seed로 하여 주변 graph를 만든다.
    
        예:
          HIF1A [nucleoplasm]
          -> component_of -> HIF-alpha, HIF complex 등
          -> located_in_pathway -> Expression of VEGFA 등
          -> participates_in_candidate_event -> pathway contained event 후보
        """
        if self.source != PathwaySource.REACTOME:
            raise ValueError("fetch_reactome_entity_neighborhood()는 source='reactome'일 때만 호출할 수 있습니다.")
    
        clean_id = self.validate_pathway_id(entity_id)
    
        if seed_payload is None:
            seed_payload = self.reactome_provider.query_reactome_id(clean_id)
    
        seed_node_id = f"REACTOME:{clean_id}"
    
        self.graph.add_node(
            seed_node_id,
            label=seed_payload.get("displayName", clean_id),
            type=seed_payload.get("schemaClass", "PhysicalEntity"),
            source="reactome",
            raw_id=clean_id,
        )
    
        added_edges = 0
        structure_count = 0
        pathway_count = 0
        event_count = 0
        reaction_detail_count = 0
    
        # 1. Entity가 포함된 Complex / Set 등 larger structures
        structures = self.reactome_provider.fetch_structures_containing_entity(clean_id)
    
        for item in structures:
            if not isinstance(item, dict):
                continue
        
            item_id = item.get("stId")
            if not item_id:
                continue
        
            item_node_id = f"REACTOME:{item_id}"
        
            self.graph.add_node(
                item_node_id,
                label=item.get("displayName", item_id),
                type=item.get("schemaClass", ""),
                source="reactome",
                raw_id=item_id,
                schemaClass=item.get("schemaClass", ""),
            )
        
            self.graph.add_edge(
                seed_node_id,
                item_node_id,
                relation_type="component_of",
                source="reactome",
            )
        
            added_edges += 1
            structure_count += 1
    
        # 2. Entity가 포함된 lower-level pathway
        pathways = self.reactome_provider.fetch_pathways_containing_entity(clean_id)
        
        for pathway in pathways:
            if not isinstance(pathway, dict):
                continue
        
            pathway_id = pathway.get("stId")
            if not pathway_id:
                continue
        
            pathway_node_id = f"REACTOME:{pathway_id}"
        
            self.graph.add_node(
                pathway_node_id,
                label=pathway.get("displayName", pathway_id),
                type=pathway.get("schemaClass", "Pathway"),
                source="reactome",
                raw_id=pathway_id,
                schemaClass=pathway.get("schemaClass", "Pathway"),
            )
        
            self.graph.add_edge(
                seed_node_id,
                pathway_node_id,
                relation_type="located_in_pathway",
                source="reactome",
            )
        
            added_edges += 1
            pathway_count += 1
    
        # 3. 관련 pathway의 contained events를 가져와 후보 event로 연결
        events = self.reactome_provider.fetch_events_containing_entity(
            clean_id,
            expand_pathway_events=True,
        )
        
        for event in events:
            if not isinstance(event, dict):
                continue
        
            event_id = event.get("stId")
            if not event_id:
                continue
        
            event_node_id = f"REACTOME:{event_id}"
            schema_class = event.get("schemaClass", "")
        
            self.graph.add_node(
                event_node_id,
                label=event.get("displayName", event_id),
                type=schema_class,
                source="reactome",
                raw_id=event_id,
                schemaClass=schema_class,
            )
        
            if not self.graph.has_edge(seed_node_id, event_node_id):
                self.graph.add_edge(
                    seed_node_id,
                    event_node_id,
                    relation_type="participates_in_candidate_event",
                    source="reactome",
                )
                added_edges += 1
        
            event_count += 1
        
            if not build_reaction_details:
                continue
        
            if schema_class not in {
                "Reaction",
                "ReactionLikeEvent",
                "BlackBoxEvent",
                "Polymerisation",
                "Depolymerisation",
            }:
                continue
        
            detail = self.reactome_provider.fetch_reaction_detail(event_id)
        
            if detail:
                result = parse_reactome_reaction_detail(detail, self.graph)
                added_edges += result.get("added_edges", 0)
                reaction_detail_count += 1
    
        return {
            "status": "success",
            "source": "reactome",
            "mode": "entity_neighborhood",
            "entity_id": clean_id,
            "schemaClass": seed_payload.get("schemaClass"),
            "displayName": seed_payload.get("displayName", clean_id),
            "structure_count": structure_count,
            "pathway_count": pathway_count,
            "event_count": event_count,
            "reaction_detail_count": reaction_detail_count,
            "added_edges": added_edges,
            "node_count": len(self.graph.nodes),
            "edge_count": len(self.graph.edges),
        }
    
    def parse_kegg_data(self, data: str) -> dict[str, Any]:
        """Compatibility wrapper for old code that passed raw KGML XML."""
        if self.source != PathwaySource.KEGG:
            raise ValueError("parse_kegg_data()는 source='kegg'일 때만 사용할 수 있습니다.")
        return parse_kegg_kgml(data, self.graph, self.graph.graph.get("pathway_id"))

    def parse_reactome_data(self, participants: list[dict[str, Any]]) -> dict[str, Any]:
        """Compatibility wrapper for old code that passed Reactome JSON lists."""
        if self.source != PathwaySource.REACTOME:
            raise ValueError("parse_reactome_data()는 source='reactome'일 때만 사용할 수 있습니다.")
        return parse_reactome_events(participants, self.graph, self.graph.graph.get("pathway_id"))

    def _parse_reaction_detail(self, reaction_id: str) -> dict[str, Any]:
        """Compatibility wrapper for the previous private method."""
        if self.source != PathwaySource.REACTOME:
            raise ValueError("_parse_reaction_detail()는 source='reactome'일 때만 사용할 수 있습니다.")
        detail = self.reactome_provider.fetch_reaction_detail(reaction_id)
        if not detail:
            return {"status": "failed", "reaction_id": reaction_id}
        return parse_reactome_reaction_detail(detail, self.graph)


    # ------------------------------------------------------------------
    # Graph utilities
    # ------------------------------------------------------------------
    def find_node_by_label(self, name_query: str) -> list[str]:
        """Find node IDs by partial label/name match."""
        query = name_query.upper().strip()
        if not query:
            return []
        return [
            node_id
            for node_id, attrs in self.graph.nodes(data=True)
            if query in str(attrs.get("label", "")).upper()
        ]

    def get_pathways_between_nodes(self, start_node, end_node, cutoff=10):
        start_node = self.normalize_node_id(start_node)
        end_node = self.normalize_node_id(end_node)
    
        if start_node not in self.graph or end_node not in self.graph:
            print(f"오류: 입력한 노드({start_node} 또는 {end_node})가 그래프에 존재하지 않습니다.")
            return []
    
        try:
            paths = list(nx.all_simple_paths(
                self.graph,
                source=start_node,
                target=end_node,
                cutoff=cutoff
            ))
            return paths
    
        except nx.NetworkXNoPath:
            return []
    
        except Exception as e:
            print(f"경로 탐색 중 오류 발생: {e}")
            return []

    def visualize_nodesedges(self, **kwargs: Any) -> None:
        """Visualize the current graph."""
        visualize_graph(self.graph, **kwargs)

    def save_graph_json(self, output_path: str | Path) -> Path:
        """Save current graph as NetworkX node-link JSON."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as file:
            json.dump(nx.node_link_data(self.graph), file, ensure_ascii=False, indent=2)
        return output_path
    
    def save_graph_summary(
        self,
        result: dict | None = None,
        output_dir: str = "results",
        node_preview_limit: int = 30,
    ) -> dict:
        """
        Graph Summary와 Node Preview를 별도 파일로 저장한다.
    
        저장 파일:
          - graph_summary.json
          - graph_summary.txt
          - node_preview.csv
        """
        from pathlib import Path
        from datetime import datetime
        import json
        import csv
    
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
        source = self.graph.graph.get("source", "")
        pathway_id = self.graph.graph.get("pathway_id", "")
        pathway_name = self.graph.graph.get("name", "")
    
        summary = {
            "source": source,
            "pathway_id": pathway_id,
            "pathway_name": pathway_name,
            "reactome_schemaClass": self.graph.graph.get("reactome_schemaClass", ""),
            "node_count": len(self.graph.nodes),
            "edge_count": len(self.graph.edges),
            "parser_status": result.get("status") if isinstance(result, dict) else None,
            "mode": result.get("mode") if isinstance(result, dict) else None,
            "result": result or {},
            "saved_at": datetime.now().isoformat(timespec="seconds"),
        }
    
        summary_json_path = output_dir / "graph_summary.json"
        summary_txt_path = output_dir / "graph_summary.txt"
        node_preview_csv_path = output_dir / "node_preview.csv"
    
        # 1. JSON 저장
        with open(summary_json_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
    
        # 2. TXT 저장
        with open(summary_txt_path, "w", encoding="utf-8") as f:
            f.write("--- Graph Summary ---\n")
            f.write(f"Source: {source}\n")
            f.write(f"Pathway ID: {pathway_id}\n")
            f.write(f"Pathway Name: {pathway_name}\n")
            f.write(f"Reactome schemaClass: {summary['reactome_schemaClass']}\n")
            f.write(f"Nodes: {summary['node_count']}\n")
            f.write(f"Edges: {summary['edge_count']}\n")
            f.write(f"Parser status: {summary['parser_status']}\n")
            f.write(f"Mode: {summary['mode']}\n")
            f.write(f"Saved at: {summary['saved_at']}\n")
    
            f.write("\n--- Node Preview ---\n")
            for node_id, attrs in list(self.graph.nodes(data=True))[:node_preview_limit]:
                f.write(
                    f"[{attrs.get('source')}] {node_id} | "
                    f"type={attrs.get('type')} | "
                    f"label={attrs.get('label')}\n"
                )
    
            remaining = len(self.graph.nodes) - node_preview_limit
            if remaining > 0:
                f.write(f"... {remaining} more nodes\n")
    
        # 3. Node Preview CSV 저장
        with open(node_preview_csv_path, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["node_id", "source", "type", "label", "raw_id"])
    
            for node_id, attrs in list(self.graph.nodes(data=True))[:node_preview_limit]:
                writer.writerow([
                    node_id,
                    attrs.get("source", ""),
                    attrs.get("type", ""),
                    attrs.get("label", ""),
                    attrs.get("raw_id", ""),
                ])
    
        return {
            "summary_json": str(summary_json_path),
            "summary_txt": str(summary_txt_path),
            "node_preview_csv": str(node_preview_csv_path),
        }

    def save_run_results(self, selected_id, result=None, paths=None, output_root="results"):
        """
        실행 결과를 results/source_pathwayid_timestamp/ 아래에 모두 저장.
        - graph JSON
        - nodes CSV
        - edges CSV
        - paths CSV
        - summary JSON
        """
        source = self.graph.graph.get(
            "source", 
            self.source.value if hasattr(self.source, "value") else self.source
        )
        source = str(source).lower()
    
        safe_id = str(selected_id).replace(":", "_").replace("/", "_")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
        run_dir = Path(output_root) / f"{source}_{safe_id}_{timestamp}"
        run_dir.mkdir(parents=True, exist_ok=True)
    
        graph_json_path = run_dir / "graph_node_link.json"
        nodes_csv_path = run_dir / "nodes.csv"
        edges_csv_path = run_dir / "edges.csv"
        paths_csv_path = run_dir / "paths.csv"
        summary_json_path = run_dir / "summary.json"
    
        # 1. graph JSON 저장
        with open(graph_json_path, "w", encoding="utf-8") as f:
            json.dump(nx.node_link_data(self.graph), f, ensure_ascii=False, indent=2)
    
        # 2. nodes.csv 저장
        with open(nodes_csv_path, "w", encoding="utf-8-sig", newline="") as f:
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
                "schemaClass",
                "compartment",
                "species",
                "reference_identifier",
                "reference_database",
                "reference_display",
            ])
        
            for node_id, attrs in self.graph.nodes(data=True):
                writer.writerow([
                    node_id,
                    attrs.get("source", source),
                    attrs.get("type", ""),
                    attrs.get("label", ""),
                    attrs.get("raw_id", ""),
                    attrs.get("compound_id", ""),
                    attrs.get("compound_name", ""),
                    attrs.get("kgml_entry_id", ""),
                    attrs.get("kgml_name", ""),
                    attrs.get("group_components", ""),
                    attrs.get("group_component_count", ""),
                    attrs.get("schemaClass", ""),
                    attrs.get("compartment", ""),
                    attrs.get("species", ""),
                    attrs.get("reference_identifier", ""),
                    attrs.get("reference_database", ""),
                    attrs.get("reference_display", ""),
                ])
    
        # 3. edges.csv 저장
        with open(edges_csv_path, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "source_node",
                "target_node",
                "relation_type",
                "source",
                "kgml_relation_type",
                "expanded_from_group",
                "original_source_node",
                "original_target_node",
            ])
        
            for u, v, attrs in self.graph.edges(data=True):
                writer.writerow([
                    u,
                    v,
                    attrs.get("relation_type", attrs.get("relation", attrs.get("rel", ""))),
                    attrs.get("source", source),
                    attrs.get("kgml_relation_type", ""),
                    attrs.get("expanded_from_group", ""),
                    attrs.get("original_source_node", ""),
                    attrs.get("original_target_node", ""),
                ])
    
        # 4. paths.csv 저장
        paths = paths or []
        with open(paths_csv_path, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["path_index", "step_index", "node_id", "label", "type"])
    
            for i, path in enumerate(paths, start=1):
                for j, node_id in enumerate(path, start=1):
                    attrs = self.graph.nodes[node_id]
                    writer.writerow([
                        i,
                        j,
                        node_id,
                        attrs.get("label", ""),
                        attrs.get("type", "")
                    ])
    
        # 5. summary 저장
        summary = {
            "source": source,
            "pathway_id": selected_id,
            "nodes": len(self.graph.nodes),
            "edges": len(self.graph.edges),
            "path_count": len(paths),
            "result": result,
            "saved_at": timestamp,
            "output_dir": str(run_dir)
        }
    
        with open(summary_json_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
    
        return run_dir

"""
    def fetch_reactome_flow(self, pathway_id: str) -> dict[str, Any]:
        Build detailed Reactome flow edges for a pathway or single reaction.
        if self.source != PathwaySource.REACTOME:
            raise ValueError("fetch_reactome_flow()는 source='reactome'일 때만 호출할 수 있습니다.")

        clean_id = self.validate_pathway_id(pathway_id)
        payload = self.reactome_provider.fetch_pathway_payload(clean_id)
        added_edges = 0
        for event in payload["events"]:
            event_id = event.get("stId")
            schema_class = event.get("schemaClass", "")
            if not event_id or schema_class not in ("ReactionLikeEvent", "Reaction"):
                continue
            detail = self.reactome_provider.fetch_reaction_detail(event_id)
            if detail:
                result = parse_reactome_reaction_detail(detail, self.graph)
                added_edges += result.get("added_edges", 0)

        return {
            "status": "success",
            "source": "reactome",
            "pathway_id": clean_id,
            "added_edges": added_edges,
            "node_count": len(self.graph.nodes),
            "edge_count": len(self.graph.edges),
        }
"""