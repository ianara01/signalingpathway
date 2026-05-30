# -*- coding: utf-8 -*-
"""Patch snippet for main.py.

Add these imports and functions to main.py. Then call run_extended_database_analysis()
after workflow_results are generated.

This patch runs:
- KEGG database coverage
- Reactome enrichment
- KOALA parsing + KO→KEGG pathway mapping
- expression mapping
- Reactome species comparison
"""

from __future__ import annotations

from pathlib import Path
import json
import subprocess
import sys


EXTENDED_ANALYSIS_CONFIG = {
    "hif1": {
        "reactome_enrichment_identifiers": [
            "HIF1A", "VEGFA", "ARNT", "VHL", "EGLN1", "EGLN2", "EGLN3",
            "EPAS1", "CREBBP", "EP300", "CA9", "EPO", "SLC2A1", "PDK1",
        ],
        "species": ["Mus musculus", "Rattus norvegicus"],
    },
    "p53": {
        "reactome_enrichment_identifiers": [
            "TP53", "MDM2", "CDKN1A", "BAX", "BBC3", "GADD45A", "ATM", "CHEK2",
        ],
        "species": ["Mus musculus", "Rattus norvegicus"],
    },
}


def _module_name(short_module: str) -> str:
    if Path("analysis").exists() and short_module.startswith("analysis."):
        return short_module
    if Path("providers").exists() and short_module.startswith("providers."):
        return short_module
    if Path("annotations").exists() and short_module.startswith("annotations."):
        return short_module
    if Path("S2SignalPathway").exists():
        return f"S2SignalPathway.{short_module}"
    return short_module


def _run_command(command: list[str], *, description: str) -> bool:
    print(f"\n--- {description} ---")
    print(" ".join(command))
    try:
        subprocess.run(command, check=True)
        return True
    except subprocess.CalledProcessError as exc:
        print(f"[warning] {description} failed: {exc}")
        return False
    except FileNotFoundError as exc:
        print(f"[warning] executable not found: {exc}")
        return False


def run_extended_database_analysis(
    target_key: str,
    workflow_results: dict,
    *,
    output_root: str = "results/analysis_extended",
    koala_file: str | None = None,
    expression_file: str | None = None,
) -> dict:
    cfg = EXTENDED_ANALYSIS_CONFIG.get(target_key, {})
    output_dir = Path(output_root) / target_key
    output_dir.mkdir(parents=True, exist_ok=True)

    kegg_dir = Path(workflow_results["kegg"]) if workflow_results.get("kegg") else None
    reactome_dir = Path(workflow_results["reactome"]) if workflow_results.get("reactome") else None

    summary = {"target_key": target_key, "output_dir": str(output_dir), "steps": {}}

    if kegg_dir and (kegg_dir / "nodes.csv").exists():
        out = output_dir / "kegg_database_coverage"
        cmd = [sys.executable, "-m", _module_name("analysis.kegg_database_coverage"), str(kegg_dir), "--output-dir", str(out)]
        ok = _run_command(cmd, description="KEGG database coverage")
        summary["steps"]["kegg_database_coverage"] = {"status": "success" if ok else "failed", "output_dir": str(out)}
    else:
        print("[skip] KEGG result directory or nodes.csv not found.")

    identifiers = cfg.get("reactome_enrichment_identifiers", [])
    if identifiers:
        out = output_dir / "reactome_enrichment"
        cmd = [sys.executable, "-m", _module_name("analysis.reactome_enrichment"), "--identifiers", *identifiers, "--output-dir", str(out), "--prefix", f"{target_key}_reactome_enrichment"]
        ok = _run_command(cmd, description="Reactome enrichment")
        summary["steps"]["reactome_enrichment"] = {"status": "success" if ok else "failed", "output_dir": str(out)}
    else:
        print("[skip] No Reactome enrichment identifiers configured.")

    if koala_file:
        annotation_out = output_dir / "annotations"
        annotation_out.mkdir(parents=True, exist_ok=True)
        parsed_ko = annotation_out / "koala_parsed.csv"
        ko_map = annotation_out / "ko_pathway_mapping.csv"
        ko_summary = annotation_out / "ko_pathway_mapping_summary.json"
        ok1 = _run_command([sys.executable, "-m", _module_name("annotations.koala_parser"), koala_file, "--output", str(parsed_ko)], description="Parse KOALA output")
        ok2 = _run_command([sys.executable, "-m", _module_name("annotations.ko_to_pathway_mapper"), str(parsed_ko), "--output-csv", str(ko_map), "--output-json", str(ko_summary)], description="Map KO to KEGG pathways")
        summary["steps"]["koala_annotation"] = {"status": "success" if ok1 and ok2 else "failed", "parsed_ko": str(parsed_ko), "ko_pathway_mapping": str(ko_map)}
    else:
        print("[skip] KOALA file not provided.")

    if expression_file:
        nodes_csv = None
        if reactome_dir and (reactome_dir / "nodes.csv").exists():
            nodes_csv = reactome_dir / "nodes.csv"
        elif kegg_dir and (kegg_dir / "nodes.csv").exists():
            nodes_csv = kegg_dir / "nodes.csv"
        if nodes_csv:
            out = output_dir / "expression_mapping"
            cmd = [sys.executable, "-m", _module_name("analysis.expression_mapping"), "--nodes", str(nodes_csv), "--expression", expression_file, "--output-dir", str(out)]
            ok = _run_command(cmd, description="Expression mapping")
            summary["steps"]["expression_mapping"] = {"status": "success" if ok else "failed", "output_dir": str(out)}
        else:
            print("[skip] nodes.csv not found for expression mapping.")
    else:
        print("[skip] expression file not provided.")

    if reactome_dir and reactome_dir.exists():
        out = output_dir / "species_comparison"
        species = cfg.get("species", ["Mus musculus", "Rattus norvegicus"])
        cmd = [sys.executable, "-m", _module_name("analysis.species_comparison"), "--reactome-result-dir", str(reactome_dir), "--species", *species, "--output-dir", str(out)]
        ok = _run_command(cmd, description="Reactome species comparison")
        summary["steps"]["species_comparison"] = {"status": "success" if ok else "failed", "output_dir": str(out)}
    else:
        print("[skip] Reactome result directory not found for species comparison.")

    summary_path = output_dir / "extended_analysis_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nExtended analysis summary saved: {summary_path}")
    return summary


# Insert this block into your existing main() after printing Workflow 완료:
"""
run_ext = input("\\nKEGG/Reactome 확장 분석을 실행하시겠습니까? (y/n): ").strip().lower()
if run_ext == "y":
    koala_file = input("KOALA 결과 파일 경로, 없으면 Enter: ").strip() or None
    expression_file = input("Expression CSV 파일 경로, 없으면 Enter: ").strip() or None
    run_extended_database_analysis(
        target_key=target_key,
        workflow_results=workflow_results,
        output_root="results/analysis_extended",
        koala_file=koala_file,
        expression_file=expression_file,
    )
"""
