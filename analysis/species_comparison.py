# -*- coding: utf-8 -*-
"""Reactome species comparison wrapper."""

from __future__ import annotations
from pathlib import Path
import json
import pandas as pd

try:
    from providers.reactome_species_provider import ReactomeSpeciesProvider
except Exception:
    from ..providers.reactome_species_provider import ReactomeSpeciesProvider


def extract_event_ids_from_reactome_result(reactome_result_dir: str | Path) -> list[str]:
    root = Path(reactome_result_dir)
    ids = []
    for name in ["reaction_event_nodes.csv", "pathway_nodes.csv", "nodes.csv"]:
        p = root / name
        if not p.exists():
            continue
        df = pd.read_csv(p, encoding="utf-8-sig")
        if "node_id" in df.columns:
            for node_id in df["node_id"].dropna().astype(str):
                if "R-HSA-" in node_id:
                    ids.append(node_id.replace("REACTOME:", ""))
    return sorted(set(ids))


def run_species_comparison(event_ids: list[str], species: list[str], output_dir: str | Path) -> dict:
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    df = ReactomeSpeciesProvider().compare_event_species(event_ids, species)
    csv_path = out / "reactome_species_comparison.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    summary = {"event_count": len(event_ids), "species": species, "output_csv": str(csv_path), "success_count": int((df["status"] == "success").sum()) if not df.empty and "status" in df.columns else 0, "failed_count": int((df["status"] == "failed").sum()) if not df.empty and "status" in df.columns else 0}
    (out / "species_comparison_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Compare Reactome event orthology across species.")
    parser.add_argument("--event-ids", nargs="*", default=None)
    parser.add_argument("--reactome-result-dir", default=None)
    parser.add_argument("--species", nargs="*", default=["Mus musculus", "Rattus norvegicus"])
    parser.add_argument("--output-dir", default="results/analysis/species_comparison")
    args = parser.parse_args()
    event_ids = args.event_ids or []
    if args.reactome_result_dir:
        event_ids.extend(extract_event_ids_from_reactome_result(args.reactome_result_dir))
    print(json.dumps(run_species_comparison(sorted(set(event_ids)), args.species, args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
