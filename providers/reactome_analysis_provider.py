# -*- coding: utf-8 -*-
"""Reactome AnalysisService provider."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import json, time, requests
import pandas as pd


class ReactomeAnalysisProvider:
    BASE_URL = "https://reactome.org/AnalysisService"

    def __init__(self, timeout: tuple[int, int] = (10, 180), sleep_sec: float = 0.2):
        self.timeout = timeout
        self.sleep_sec = sleep_sec

    @staticmethod
    def as_text_payload(identifiers: list[str], header: str = "#Genes") -> str:
        return "\n".join([header] + [str(x).strip() for x in identifiers if str(x).strip()]) + "\n"

    def submit_identifiers_projection(
        self,
        identifiers: list[str] | str,
        *,
        page_size: int = 100,
        page: int = 1,
        species: str = "Homo sapiens",
        include_disease: bool = True,
        p_value: float = 1.0,
        resource: str = "TOTAL",
    ) -> dict[str, Any]:
        payload = identifiers if isinstance(identifiers, str) else self.as_text_payload(identifiers)
        url = f"{self.BASE_URL}/identifiers/projection"
        params = {
            "pageSize": page_size,
            "page": page,
            "species": species,
            "includeDisease": str(include_disease).lower(),
            "pValue": p_value,
            "resource": resource,
        }
        response = requests.post(
            url, params=params, data=payload.encode("utf-8"), timeout=self.timeout,
            headers={"Content-Type": "text/plain", "Accept": "application/json", "User-Agent": "S2SignalPathway/1.0"},
        )
        response.raise_for_status()
        time.sleep(self.sleep_sec)
        return response.json()

    def submit_file_projection(self, input_file: str | Path, **kwargs: Any) -> dict[str, Any]:
        return self.submit_identifiers_projection(Path(input_file).read_text(encoding="utf-8"), **kwargs)

    def token_result(self, token: str, *, page_size: int = 100, page: int = 1, resource: str = "TOTAL", species: str = "Homo sapiens") -> dict[str, Any]:
        response = requests.get(
            f"{self.BASE_URL}/token/{token}",
            params={"pageSize": page_size, "page": page, "resource": resource, "species": species},
            timeout=self.timeout,
            headers={"Accept": "application/json", "User-Agent": "S2SignalPathway/1.0"},
        )
        response.raise_for_status()
        return response.json()

    @staticmethod
    def flatten_pathway_results(result: dict[str, Any]) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []

        def walk(obj: Any, context: dict[str, Any] | None = None) -> None:
            context = context or {}
            if isinstance(obj, dict):
                st_id = obj.get("stId") or obj.get("stIdVersion") or obj.get("identifier")
                name = obj.get("name") or obj.get("displayName")
                entities = obj.get("entities") or {}
                reactions = obj.get("reactions") or {}
                if st_id or name:
                    rows.append({
                        "pathway_id": st_id or "",
                        "pathway_name": name or "",
                        "species": obj.get("speciesName") or obj.get("species") or context.get("species", ""),
                        "entities_found": entities.get("found", "") if isinstance(entities, dict) else "",
                        "entities_total": entities.get("total", "") if isinstance(entities, dict) else "",
                        "entities_ratio": entities.get("ratio", "") if isinstance(entities, dict) else "",
                        "reactions_found": reactions.get("found", "") if isinstance(reactions, dict) else "",
                        "reactions_total": reactions.get("total", "") if isinstance(reactions, dict) else "",
                        "p_value": obj.get("pValue", ""),
                        "fdr": obj.get("fdr", ""),
                        "resource": context.get("resource", ""),
                    })
                for key, value in obj.items():
                    ctx = dict(context)
                    if key.lower() == "resource":
                        ctx["resource"] = str(value)
                    walk(value, ctx)
            elif isinstance(obj, list):
                for item in obj:
                    walk(item, context)

        walk(result)
        if not rows:
            return pd.DataFrame()
        df = pd.DataFrame(rows).drop_duplicates()
        if "fdr" in df.columns:
            df["_fdr"] = pd.to_numeric(df["fdr"], errors="coerce")
            df = df.sort_values(["_fdr", "pathway_name"], na_position="last").drop(columns=["_fdr"])
        return df.reset_index(drop=True)

    def run_identifier_enrichment_to_files(self, identifiers: list[str], output_dir: str | Path, prefix: str = "reactome_enrichment", **kwargs: Any) -> dict[str, Any]:
        out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
        result = self.submit_identifiers_projection(identifiers, **kwargs)
        json_path, csv_path = out / f"{prefix}.json", out / f"{prefix}.csv"
        json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        df = self.flatten_pathway_results(result)
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")
        return {"json": str(json_path), "csv": str(csv_path), "row_count": int(len(df))}


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Run Reactome enrichment/projection.")
    parser.add_argument("--identifiers", nargs="*", default=None)
    parser.add_argument("--input-file", default=None)
    parser.add_argument("--output-dir", default="results/reactome_analysis")
    parser.add_argument("--prefix", default="reactome_enrichment")
    args = parser.parse_args()
    provider = ReactomeAnalysisProvider()
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    result = provider.submit_file_projection(args.input_file) if args.input_file else provider.submit_identifiers_projection(args.identifiers or [])
    (out / f"{args.prefix}.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    provider.flatten_pathway_results(result).to_csv(out / f"{args.prefix}.csv", index=False, encoding="utf-8-sig")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
