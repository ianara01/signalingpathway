# -*- coding: utf-8 -*-
"""Map KO identifiers to KEGG pathways using KEGG REST link/pathway/ko."""

from __future__ import annotations
from pathlib import Path
import json
import requests
import pandas as pd

KEGG_REST = "https://rest.kegg.jp"


def map_ko_to_pathways(ko_ids: list[str], timeout: int = 30) -> pd.DataFrame:
    rows = []
    for ko in sorted(set([str(x).replace("ko:", "").strip() for x in ko_ids if str(x).strip()])):
        try:
            r = requests.get(f"{KEGG_REST}/link/pathway/ko:{ko}", timeout=timeout, headers={"User-Agent": "S2SignalPathway/1.0"})
            r.raise_for_status()
            for line in r.text.splitlines():
                parts = line.split("\t")
                if len(parts) == 2:
                    rows.append({"ko_id": ko, "pathway_id": parts[1].replace("path:", ""), "source": parts[0]})
        except Exception as exc:
            rows.append({"ko_id": ko, "pathway_id": "", "source": f"ko:{ko}", "error": str(exc)})
    return pd.DataFrame(rows)


def map_ko_file_to_pathways(input_file: str | Path, output_csv: str | Path | None = None) -> pd.DataFrame:
    path = Path(input_file)
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path, encoding="utf-8-sig")
        if "ko_id" not in df.columns:
            raise ValueError("CSV must contain ko_id column.")
        ko_ids = df["ko_id"].dropna().astype(str).tolist()
    else:
        ko_ids = [x.strip() for x in path.read_text(encoding="utf-8", errors="ignore").splitlines() if x.strip()]
    mapped = map_ko_to_pathways(ko_ids)
    if output_csv:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        mapped.to_csv(output_csv, index=False, encoding="utf-8-sig")
    return mapped


def summarize_ko_pathway_mapping(mapped: pd.DataFrame) -> dict:
    return {
        "ko_count": int(mapped["ko_id"].nunique()) if not mapped.empty and "ko_id" in mapped.columns else 0,
        "pathway_count": int(mapped["pathway_id"].nunique()) if not mapped.empty and "pathway_id" in mapped.columns else 0,
        "top_pathways": mapped["pathway_id"].value_counts().head(20).to_dict() if not mapped.empty and "pathway_id" in mapped.columns else {},
    }


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Map KO IDs to KEGG pathways.")
    parser.add_argument("input_file")
    parser.add_argument("--output-csv", default="results/annotations/ko_pathway_mapping.csv")
    parser.add_argument("--output-json", default="results/annotations/ko_pathway_mapping_summary.json")
    args = parser.parse_args()
    mapped = map_ko_file_to_pathways(args.input_file, args.output_csv)
    summary = summarize_ko_pathway_mapping(mapped)
    Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_json).write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Saved: {args.output_csv}")
    print(f"Saved: {args.output_json}")


if __name__ == "__main__":
    main()
