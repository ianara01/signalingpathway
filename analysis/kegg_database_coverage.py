# -*- coding: utf-8 -*-
"""KEGG database coverage analysis."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import json
import pandas as pd

try:
    from providers.kegg_extended_provider import KeggExtendedProvider
except Exception:
    from ..providers.kegg_extended_provider import KeggExtendedProvider


def _read_csv(path: str | Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, encoding="utf-8-sig")
    except UnicodeDecodeError:
        return pd.read_csv(path)


def extract_kegg_genes_and_compounds(nodes: pd.DataFrame) -> tuple[list[str], list[str]]:
    genes, compounds = [], []
    if nodes.empty:
        return genes, compounds
    for _, row in nodes.iterrows():
        ntype = str(row.get("type", "")).lower()
        raw_id = str(row.get("raw_id", "")).strip()
        node_id = str(row.get("node_id", "")).replace("KEGG:", "").strip()
        value = raw_id or node_id
        if ntype == "gene" and value:
            if value.isdigit() or value.startswith("hsa:"):
                genes.append(value)
        elif ntype == "compound" and value:
            compounds.append(value.replace("cpd:", ""))
    return sorted(set(genes)), sorted(set(compounds))


def run_kegg_database_coverage(kegg_result_dir: str | Path, output_dir: str | Path, *, max_genes: int | None = None, max_compounds: int | None = None) -> dict[str, Any]:
    result_dir, out = Path(kegg_result_dir), Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    nodes = _read_csv(result_dir / "nodes.csv")
    genes, compounds = extract_kegg_genes_and_compounds(nodes)
    if max_genes:
        genes = genes[:max_genes]
    if max_compounds:
        compounds = compounds[:max_compounds]

    provider = KeggExtendedProvider()
    summary: dict[str, Any] = {"kegg_result_dir": str(result_dir), "output_dir": str(out), "gene_count": len(genes), "compound_count": len(compounds)}

    if genes:
        cov = provider.build_gene_coverage(genes)
        cov["gene_annotations"].to_csv(out / "kegg_gene_annotations.csv", index=False, encoding="utf-8-sig")
        cov["links"].to_csv(out / "kegg_gene_database_links.csv", index=False, encoding="utf-8-sig")
        summary["gene_link_target_database_counts"] = cov["links"]["target_database"].value_counts().to_dict() if not cov["links"].empty and "target_database" in cov["links"].columns else {}
    if compounds:
        cov = provider.build_compound_coverage(compounds)
        cov["compound_annotations"].to_csv(out / "kegg_compound_annotations.csv", index=False, encoding="utf-8-sig")
        cov["links"].to_csv(out / "kegg_compound_database_links.csv", index=False, encoding="utf-8-sig")
        summary["compound_link_target_database_counts"] = cov["links"]["target_database"].value_counts().to_dict() if not cov["links"].empty and "target_database" in cov["links"].columns else {}

    coverage_matrix = [
        {"category": "KEGG PATHWAY", "covered": True, "evidence_file": "nodes.csv/edges.csv"},
        {"category": "KEGG GENES", "covered": bool(genes), "evidence_file": "kegg_gene_annotations.csv"},
        {"category": "KEGG COMPOUND", "covered": bool(compounds), "evidence_file": "kegg_compound_annotations.csv"},
        {"category": "KEGG DRUG", "covered": bool(summary.get("gene_link_target_database_counts", {}).get("drug", 0)), "evidence_file": "kegg_gene_database_links.csv"},
        {"category": "KEGG DISEASE", "covered": bool(summary.get("gene_link_target_database_counts", {}).get("disease", 0)), "evidence_file": "kegg_gene_database_links.csv"},
        {"category": "KEGG ENZYME", "covered": bool(summary.get("gene_link_target_database_counts", {}).get("enzyme", 0) or summary.get("compound_link_target_database_counts", {}).get("enzyme", 0)), "evidence_file": "kegg_*_database_links.csv"},
        {"category": "KEGG REACTION", "covered": bool(summary.get("gene_link_target_database_counts", {}).get("reaction", 0) or summary.get("compound_link_target_database_counts", {}).get("reaction", 0)), "evidence_file": "kegg_*_database_links.csv"},
    ]
    pd.DataFrame(coverage_matrix).to_csv(out / "kegg_database_coverage_matrix.csv", index=False, encoding="utf-8-sig")
    summary["coverage_matrix"] = coverage_matrix
    (out / "kegg_database_coverage_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Analyze KEGG database coverage from KEGG result folder.")
    parser.add_argument("kegg_result_dir")
    parser.add_argument("--output-dir", default="results/analysis/kegg_database_coverage")
    parser.add_argument("--max-genes", type=int, default=None)
    parser.add_argument("--max-compounds", type=int, default=None)
    args = parser.parse_args()
    print(json.dumps(run_kegg_database_coverage(args.kegg_result_dir, args.output_dir, max_genes=args.max_genes, max_compounds=args.max_compounds), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
