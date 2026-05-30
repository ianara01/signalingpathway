# -*- coding: utf-8 -*-
"""Expression mapping onto parsed pathway nodes."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import json
import pandas as pd


def _read_csv(path: str | Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, encoding="utf-8-sig")
    except UnicodeDecodeError:
        return pd.read_csv(path)


def _first_existing(columns: list[str], candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def map_expression_to_nodes(nodes: pd.DataFrame, expression: pd.DataFrame, *, expr_id_col: str | None = None) -> pd.DataFrame:
    if nodes.empty or expression.empty:
        return pd.DataFrame()
    expr_id_col = expr_id_col or _first_existing(list(expression.columns), ["gene", "symbol", "identifier", "id", "target", "gene_symbol"]) or expression.columns[0]
    logfc_col = _first_existing(list(expression.columns), ["log2FC", "logFC", "fold_change", "fc"])
    pval_col = _first_existing(list(expression.columns), ["pvalue", "p_value", "pval", "padj", "fdr"])
    expr = expression.copy()
    expr["_match_key"] = expr[expr_id_col].fillna("").astype(str).str.upper()
    search_cols = [c for c in ["label", "raw_id", "node_id", "reference_identifier", "reference_display"] if c in nodes.columns]
    rows = []
    for _, node in nodes.iterrows():
        haystack = " ".join(str(node.get(c, "")) for c in search_cols).upper()
        matches = expr[expr["_match_key"].apply(lambda x: bool(x and x in haystack))]
        for _, m in matches.iterrows():
            row = node.to_dict()
            row.update({"expression_identifier": m.get(expr_id_col, ""), "log2FC": m.get(logfc_col, "") if logfc_col else "", "p_value_or_fdr": m.get(pval_col, "") if pval_col else ""})
            rows.append(row)
    return pd.DataFrame(rows)


def run_expression_mapping(nodes_csv: str | Path, expression_csv: str | Path, output_dir: str | Path) -> dict[str, Any]:
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    overlay = map_expression_to_nodes(_read_csv(nodes_csv), _read_csv(expression_csv))
    overlay_path = out / "node_expression_overlay.csv"
    overlay.to_csv(overlay_path, index=False, encoding="utf-8-sig")
    summary = {"nodes_csv": str(nodes_csv), "expression_csv": str(expression_csv), "overlay_rows": int(len(overlay)), "matched_node_count": int(overlay["node_id"].nunique()) if not overlay.empty and "node_id" in overlay.columns else 0, "output_csv": str(overlay_path)}
    (out / "expression_mapping_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Map expression CSV onto pathway nodes.csv.")
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--output-dir", default="results/analysis/expression_mapping")
    args = parser.parse_args()
    print(json.dumps(run_expression_mapping(args.nodes, args.expression, args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
