# -*- coding: utf-8 -*-
"""Reactome species/orthology provider."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import json, time, requests
import pandas as pd


class ReactomeSpeciesProvider:
    BASE_URL = "https://reactome.org/ContentService"

    def __init__(self, timeout: tuple[int, int] = (10, 120), sleep_sec: float = 0.2):
        self.timeout = timeout
        self.sleep_sec = sleep_sec

    def _get_json(self, endpoint: str, params: dict[str, Any] | None = None) -> Any:
        url = endpoint if endpoint.startswith("http") else f"{self.BASE_URL}{endpoint}"
        response = requests.get(url, params=params or {}, timeout=self.timeout, headers={"Accept": "application/json", "User-Agent": "S2SignalPathway/1.0"})
        response.raise_for_status()
        time.sleep(self.sleep_sec)
        return response.json()

    def event_hierarchy_for_species(self, species: str = "9606") -> Any:
        return self._get_json(f"/data/eventsHierarchy/{species}")

    def orthology_for_event(self, reactome_id: str, species: str) -> dict[str, Any]:
        candidates = [
            f"/data/orthology/{reactome_id}/species/{species}",
            f"/data/orthology/{reactome_id}/{species}",
            f"/data/event/{reactome_id}/orthologous/{species}",
        ]
        errors = []
        for endpoint in candidates:
            try:
                return {"reactome_id": reactome_id, "species": species, "endpoint": endpoint, "status": "success", "data": self._get_json(endpoint)}
            except Exception as exc:
                errors.append(f"{endpoint}: {exc}")
        return {"reactome_id": reactome_id, "species": species, "status": "failed", "errors": errors}

    @staticmethod
    def flatten_event_hierarchy(hierarchy: Any) -> pd.DataFrame:
        rows = []
        def walk(node: Any, parent_id: str = "", depth: int = 0) -> None:
            if isinstance(node, dict):
                st_id = node.get("stId") or node.get("identifier") or node.get("stableIdentifier", {}).get("identifier", "")
                rows.append({
                    "event_id": st_id,
                    "event_name": node.get("name") or node.get("displayName", ""),
                    "type": node.get("type", "") or node.get("schemaClass", ""),
                    "species": node.get("species", ""),
                    "parent_id": parent_id,
                    "depth": depth,
                })
                for child in node.get("children") or node.get("hasEvent") or node.get("events") or []:
                    walk(child, st_id, depth + 1)
            elif isinstance(node, list):
                for item in node:
                    walk(item, parent_id, depth)
        walk(hierarchy)
        return pd.DataFrame(rows)

    def compare_event_species(self, event_ids: list[str], species_list: list[str]) -> pd.DataFrame:
        rows = []
        for event_id in event_ids:
            for species in species_list:
                result = self.orthology_for_event(event_id, species)
                data = result.get("data")
                count = len(data) if isinstance(data, list) else (1 if isinstance(data, dict) else 0)
                rows.append({
                    "event_id": event_id,
                    "species": species,
                    "status": result.get("status"),
                    "ortholog_count": count,
                    "endpoint": result.get("endpoint", ""),
                    "error": "; ".join(result.get("errors", [])) if result.get("errors") else "",
                })
        return pd.DataFrame(rows)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Reactome species comparison helper.")
    parser.add_argument("--event-ids", nargs="*", default=[])
    parser.add_argument("--species", nargs="*", default=["Mus musculus", "Rattus norvegicus"])
    parser.add_argument("--output-dir", default="results/reactome_species")
    args = parser.parse_args()
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    df = ReactomeSpeciesProvider().compare_event_species(args.event_ids, args.species)
    df.to_csv(out / "species_comparison.csv", index=False, encoding="utf-8-sig")
    (out / "species_summary.json").write_text(json.dumps({"event_ids": args.event_ids, "species": args.species}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
