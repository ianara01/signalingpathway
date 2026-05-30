# -*- coding: utf-8 -*-
"""Reactome enrichment analysis wrapper."""

from __future__ import annotations
from pathlib import Path
import json
import pandas as pd

try:
    from providers.reactome_analysis_provider import ReactomeAnalysisProvider
except Exception:
    from ..providers.reactome_analysis_provider import ReactomeAnalysisProvider


def load_identifiers_from_file(path: str | Path) -> list[str]:
    p = Path(path)
    if p.suffix.lower() == ".csv":
        df = pd.read_csv(p, encoding="utf-8-sig")
        for col in ["gene", "symbol", "identifier", "id", "target", "query_id"]:
            if col in df.columns:
                return df[col].dropna().astype(str).tolist()
        return df.iloc[:, 0].dropna().astype(str).tolist()
    return [line.strip() for line in p.read_text(encoding="utf-8", errors="ignore").splitlines() if line.strip() and not line.startswith("#")]


def run_reactome_enrichment(identifiers: list[str], output_dir: str | Path, *, prefix: str = "reactome_enrichment") -> dict:
    return ReactomeAnalysisProvider().run_identifier_enrichment_to_files(identifiers, output_dir, prefix=prefix)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Run Reactome pathway enrichment.")
    parser.add_argument("--identifiers", nargs="*", default=None)
    parser.add_argument("--input-file", default=None)
    parser.add_argument("--output-dir", default="results/analysis/reactome_enrichment")
    parser.add_argument("--prefix", default="reactome_enrichment")
    args = parser.parse_args()
    ids = args.identifiers or []
    if args.input_file:
        ids.extend(load_identifiers_from_file(args.input_file))
    print(json.dumps(run_reactome_enrichment(ids, args.output_dir, prefix=args.prefix), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
