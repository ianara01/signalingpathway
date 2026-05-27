"""
Reactome ContentService provider.

This module is intentionally Reactome-only. It returns raw JSON payloads and
keeps graph parsing in ``parsers/reactome_json_parser.py``.

Created on Mon Sep  8 14:15:50 2025

@author: USER, SANG JIN PARK
"""

from __future__ import annotations

from typing import Any

from pathlib import Path
import requests
import time

class ReactomeProvider:
    """Thin client for Reactome ContentService API."""

    BASE_URL = "https://reactome.org/ContentService/data"

    def __init__(self, timeout: int | tuple[int, int] = (10, 120)):
        self.timeout = timeout
    
    def _as_list(self, data: Any) -> list[Any]:
        """None/dict/list/scalar 응답을 list로 정규화."""
        if data is None:
            return []
        if isinstance(data, list):
            return data
        return [data]
    
    
    def _is_stable_id(self, value: Any) -> bool:
        return isinstance(value, str) and value.upper().startswith("R-")
    
    
    def _safe_fetch_if_stable_id(self, value: Any) -> dict[str, Any] | None:
        """
        응답 항목이 stable ID 문자열이면 /query/{id}로 확장.
        int DB_ID는 stable ID가 아니므로 여기서는 무시.
        """
        if not self._is_stable_id(value):
            return None
    
        try:
            return self.query_reactome_id(value)
        except Exception:
            return None
    
    
    def normalize_reactome_items(self, data: Any) -> list[dict[str, Any]]:
        """Normalize Reactome response into list[dict].

        dict items are kept, stable-id strings are expanded via /query/{id},
        and int DB_ID / None / other scalar items are ignored to prevent
        `'int' object has no attribute 'get'` errors.
        
        Reactome API 응답을 list[dict]로 정규화.
    
        - dict 중 stId가 있으면 그대로 사용
        - list는 내부 항목별 처리
        - str stable ID, 예: R-HSA-xxxx 는 query로 확장
        - int DB_ID는 .get() 에러 방지를 위해 제외
        """
        normalized: list[dict[str, Any]] = []
    
        for item in self._as_list(data):
            if isinstance(item, dict):
                if item.get("stId") or item.get("dbId"):
                    normalized.append(item)
                continue
    
            expanded = self._safe_fetch_if_stable_id(item)
            if expanded:
                normalized.append(expanded)
    
            # int, float, None 등은 여기서 무시
    
        return normalized

    @staticmethod
    def normalize_pathway_id(pathway_id: str) -> str:
        """Normalize Reactome stable IDs such as ``R-HSA-168256``.
        """
        if not pathway_id:
            raise ValueError("Reactome pathway_id is empty.")

        clean_id = str(pathway_id).strip().upper()

        if clean_id.startswith("REACTOME:"):
            clean_id = clean_id.replace("REACTOME:", "", 1)

        if not clean_id.startswith("R-"):
            raise ValueError(
                "Reactome pathway/reaction/entity ID 형식이 아닙니다. "
                "예: R-HSA-168256"
            )
            
        return clean_id

    def _get_json(
        self,
        endpoint: str,
        *,
        default: Any = None,
        raise_on_error: bool = False,
        retries: int = 3,
        backoff: int = 5,
    ) -> Any:
        """
        Reactome ContentService GET helper.

        endpoint 예:
          /query/R-HSA-1234135
          /pathway/R-HSA-1234174/containedEvents
        """
        
        url = f"{self.BASE_URL}{endpoint}"
        last_error = None
    
        for attempt in range(1, retries + 1):
            try:
                response = requests.get(
                    url,
                    timeout=self.timeout,
                    headers={
                        "accept": "application/json",
                        "User-Agent": "S2SignalPathway/1.0",
                    },
                )
                response.raise_for_status()
                return response.json()
    
            except requests.exceptions.ConnectTimeout:
                last_error = TimeoutError(
                    f"Reactome 서버 연결 시간이 초과되었습니다. "
                    f"timeout={self.timeout}, url={url}"
                )
    
            except requests.exceptions.ReadTimeout:
                last_error = TimeoutError(
                    f"Reactome 서버 응답 대기 시간이 초과되었습니다. "
                    f"timeout={self.timeout}, url={url}"
                )
    
            except requests.RequestException as exc:
                last_error = ConnectionError(
                    f"Reactome 요청 중 네트워크 오류가 발생했습니다: {url} ({exc})"
                )
    
            time.sleep(backoff * attempt)
    
        if raise_on_error:
            raise last_error
    
        return default


    # ---------------------------------------------------------------------
    # Basic object query
    # ---------------------------------------------------------------------

    def query_reactome_id(self, stable_id: str) -> dict[str, Any]:
        """
        Fetch a single Reactome object by stable ID.

        Reactome ID가 Pathway인지, ReactionLikeEvent인지, PhysicalEntity인지
        schemaClass를 확인하기 위한 기본 진입점.
        """
        clean_id = self.normalize_pathway_id(stable_id)
        data = self._get_json(
            f"/query/{clean_id}",
            default=None,
            raise_on_error=True,
        )

        if not isinstance(data, dict):
            raise ValueError(f"Unexpected Reactome response for {clean_id}")

        return data

    def fetch_query(self, stable_id: str) -> dict[str, Any]:
        """
        Backward-compatible alias.

        기존 코드가 fetch_query()를 사용하고 있다면 그대로 동작하게 유지.
        """
        return self.query_reactome_id(stable_id)

    # ---------------------------------------------------------------------
    # Pathway / Event queries
    # ---------------------------------------------------------------------

    def fetch_contained_events(self, pathway_id: str) -> list[dict[str, Any]] | None:
        clean_id = self.normalize_pathway_id(pathway_id)
        data = self._get_json(
            f"/pathway/{clean_id}/containedEvents",
            default=None,
            raise_on_error=False,
        )
    
        if data is None:
            return None
    
        return self.normalize_reactome_items(data)

    def fetch_pathway_payload(self, pathway_id: str) -> dict[str, Any]:
        """
        Return pathway payload.

        주의:
          이 함수는 Pathway 중심 payload입니다.
          Entity / Reaction 자동 분기는 SignalingPathway.fetch_reactome()에서 합니다.
        """
        clean_id = self.normalize_pathway_id(pathway_id)
        query = self.query_reactome_id(clean_id)
        events = self.fetch_contained_events(clean_id)

        if events is None:
            events = []

        return {
            "source": "reactome",
            "pathway_id": clean_id,
            "format": "json",
            "schemaClass": query.get("schemaClass"),
            "displayName": query.get("displayName"),
            "query": query,
            "events": events,
        }

    def fetch_reaction_detail(self, reaction_id: str) -> dict[str, Any] | None:
        """
        Fetch detailed ReactionLikeEvent information.

        실패 시 None 반환.
        """
        try:
            return self.query_reactome_id(reaction_id)
        except Exception:
            return None

    def fetch_event_ancestors(self, event_id: str) -> list[Any]:
        """
        Fetch all ancestor paths for an event.

        Reaction/Event가 어떤 pathway hierarchy에 포함되는지 확인할 때 사용.
        """
        clean_id = self.normalize_pathway_id(event_id)
        data = self._get_json(
            f"/event/{clean_id}/ancestors",
            default=[],
            raise_on_error=False,
        )
        return data if isinstance(data, list) else [data]

    # ---------------------------------------------------------------------
    # PhysicalEntity queries
    # ---------------------------------------------------------------------

    def fetch_structures_containing_entity(
        self,
        entity_id: str,
    ) -> list[dict[str, Any]]:
        """
        PhysicalEntity가 component로 포함된 Complex / Set 등 larger structures 조회.
        응답에 int DB_ID가 섞여도 list[dict]만 반환하도록 정규화.
        """
        clean_id = self.normalize_pathway_id(entity_id)
    
        data = self._get_json(
            f"/entity/{clean_id}/componentOf",
            default=[],
            raise_on_error=False,
        )
    
        return self.normalize_reactome_items(data)
    
    
    def fetch_other_forms_of_entity(
        self,
        entity_id: str,
    ) -> list[dict[str, Any]]:
        clean_id = self.normalize_pathway_id(entity_id)
    
        data = self._get_json(
            f"/entity/{clean_id}/otherForms",
            default=[],
            raise_on_error=False,
        )
    
        return self.normalize_reactome_items(data)
    
    
    def fetch_pathways_containing_entity(
        self,
        entity_id: str,
        *,
        diagram_only: bool = False,
        all_forms: bool = False,
    ) -> list[dict[str, Any]]:
        """
        PhysicalEntity가 포함된 lower-level pathway 목록 조회.
        응답을 반드시 list[dict]로 정규화.
        """
        clean_id = self.normalize_pathway_id(entity_id)
    
        if diagram_only and all_forms:
            endpoint = f"/pathways/low/diagram/entity/{clean_id}/allForms"
        elif diagram_only:
            endpoint = f"/pathways/low/diagram/entity/{clean_id}"
        else:
            endpoint = f"/pathways/low/entity/{clean_id}"
    
        data = self._get_json(
            endpoint,
            default=[],
            raise_on_error=False,
        )
    
        return self.normalize_reactome_items(data)
    
    
    def fetch_events_containing_entity(
        self,
        entity_id: str,
        *,
        expand_pathway_events: bool = True,
    ) -> list[dict[str, Any]]:
        """
        Entity가 포함된 lower-level pathways를 가져온 뒤,
        각 pathway의 containedEvents를 가져온다.
        모든 단계에서 dict 항목만 유지한다.
        """
        clean_id = self.normalize_pathway_id(entity_id)
    
        pathways = self.fetch_pathways_containing_entity(clean_id)
    
        events: list[dict[str, Any]] = []
        seen: set[str] = set()
    
        for pathway in pathways:
            if not isinstance(pathway, dict):
                continue
    
            pathway_id = pathway.get("stId")
            if not pathway_id:
                continue
    
            if pathway_id not in seen:
                events.append(pathway)
                seen.add(pathway_id)
    
            if not expand_pathway_events:
                continue
    
            contained = self.fetch_contained_events(pathway_id)
            contained = self.normalize_reactome_items(contained)
    
            for event in contained:
                event_id = event.get("stId")
                if event_id and event_id not in seen:
                    events.append(event)
                    seen.add(event_id)
    
        return events


    # ---------------------------------------------------------------------
    # Flatten helpers
    # ---------------------------------------------------------------------

    def flatten_reactome_tree(self, data: Any) -> list[dict[str, Any]]:
        """
        Reactome tree/list/dict 응답을 stId 기준 list로 flatten.
        """
        flattened: list[dict[str, Any]] = []

        def walk(obj: Any) -> None:
            if isinstance(obj, dict):
                if obj.get("stId"):
                    flattened.append(obj)

                for key in (
                    "children",
                    "hasEvent",
                    "events",
                    "pathway",
                    "ancestors",
                    "parents",
                ):
                    value = obj.get(key)
                    if value:
                        walk(value)

            elif isinstance(obj, list):
                for item in obj:
                    walk(item)

        walk(data)
        return flattened

    # ---------------------------------------------------------------------
    # recursive expansion
    # ---------------------------------------------------------------------
    def fetch_contained_events_recursive(
        self,
        pathway_id: str,
        *,
        max_depth: int = 3,
    ) -> list[dict[str, Any]]:
        clean_id = self.normalize_pathway_id(pathway_id)
    
        collected: list[dict[str, Any]] = []
        seen: set[str] = set()
    
        def walk(current_id: str, depth: int) -> None:
            if depth > max_depth:
                return
    
            events = self.fetch_contained_events(current_id) or []
            events = self.normalize_reactome_items(events)
    
            for event in events:
                event_id = event.get("stId")
                if not event_id or event_id in seen:
                    continue
    
                seen.add(event_id)
                collected.append(event)
    
                schema_class = event.get("schemaClass", "")
                if schema_class in {"Pathway", "TopLevelPathway"}:
                    walk(event_id, depth + 1)
    
        walk(clean_id, 1)
        return collected
    
    # ---------------------------------------------------------------------
    # SBML extraction
    # ---------------------------------------------------------------------
    def export_event_sbml(
        self,
        reactome_id: str,
        output_path: str | Path,
        *,
        overwrite: bool = False,
    ) -> Path:
        """
        Reactome event/pathway/reaction ID를 SBML로 export하여 저장한다.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
        if output_path.exists() and not overwrite:
            return output_path
    
        url = f"{self.BASE_URL}/exporter/event/{reactome_id}.sbml"
    
        response = requests.get(
            url,
            timeout=self.timeout,
            headers={
                "Accept": "application/xml, text/xml, application/sbml+xml, */*",
                "User-Agent": "S2SignalPathway/1.0",
            },
        )
        response.raise_for_status()
    
        text = response.text
    
        if "<sbml" not in text[:2000].lower():
            raise ValueError(
                f"Downloaded content does not look like SBML: {reactome_id}"
            )
    
        output_path.write_text(text, encoding="utf-8")
        return output_path
    