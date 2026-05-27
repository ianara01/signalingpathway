# -*- coding: utf-8 -*-
"""Model evidence checker.

This module inspects parsed S2SignalPathway outputs and determines what kind of
mathematical modeling evidence exists for a target.

Evidence levels:
- topology_graph_model: nodes and directed labeled edges exist
- signed_regulatory_model_candidate: activation/inhibition/expression edges exist
- boolean_model_candidate: signed regulatory inputs exist for at least one node
- reactome_reaction_event_model: reaction/event nodes and input/output edges exist
- stoichiometric_like_model_candidate: input/output reaction table includes chemical entities
- sbml_candidate: Reactome reaction model exists or SBML files are present
- kinetic_ode_model_evidence: SBML with kineticLaw/parameters/rules exists
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import json
import pandas as pd

try:
    import networkx as nx
except Exception:  # pragma: no cover
    nx = None

from .target_type_classifier import classify_targets_by_keyword, classify_nodes_dataframe


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, encoding="utf-8-sig")
    except UnicodeDecodeError:
        return pd.read_csv(path)


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    return data if isinstance(data, dict) else {}


def _find_first(result_dir: Path, *names: str) -> pd.DataFrame:
    for name in names:
        p = result_dir / name
        if p.exists():
            return _read_csv(p)
    return pd.DataFrame()


def load_result_tables(result_dir: str | Path) -> dict[str, Any]:
    """Load common S2SignalPathway result CSVs from one result directory."""
    root = Path(result_dir)
    csv_tables: dict[str, pd.DataFrame] = {}
    for csv_path in sorted(root.glob("*.csv")):
        csv_tables[csv_path.name] = _read_csv(csv_path)

    return {
        "root": str(root),
        "csv_tables": csv_tables,
        "nodes": csv_tables.get("nodes.csv", pd.DataFrame()),
        "edges": csv_tables.get("edges.csv", pd.DataFrame()),
        "reaction_io": csv_tables.get("reaction_io.csv", pd.DataFrame()),
        "event_order_links": csv_tables.get("event_order_links.csv", pd.DataFrame()),
        "pathway_event_links": csv_tables.get("pathway_event_links.csv", pd.DataFrame()),
        "chemical_nodes": csv_tables.get("chemical_nodes.csv", pd.DataFrame()),
        "complex_nodes": csv_tables.get("complex_nodes.csv", pd.DataFrame()),
        "entity_pathway_links": csv_tables.get("entity_pathway_links.csv", pd.DataFrame()),
        "summary": _read_json(root / "workflow_summary.json"),
        "graph_json": _read_json(root / "graph_node_link.json"),
    }


def _contains(df: pd.DataFrame, column: str, keyword: str) -> pd.Series:
    if df.empty or column not in df.columns:
        return pd.Series(False, index=df.index)
    return df[column].fillna("").astype(str).str.contains(keyword, case=False, regex=False)


def _value_counts(df: pd.DataFrame, column: str) -> dict[str, int]:
    if df.empty or column not in df.columns:
        return {}
    return {str(k): int(v) for k, v in df[column].fillna("").value_counts().items()}


def _relation_polarity(relation_type: str) -> int | None:
    r = str(relation_type).lower()
    if "inhibition" in r or "repression" in r or "negative" in r or "ubiquitination" in r:
        return -1
    if "activation" in r or "expression" in r or "positive" in r or "phosphorylation" in r:
        return 1
    return None


def _signed_edge_stats(edges: pd.DataFrame) -> dict[str, Any]:
    if edges.empty or "relation_type" not in edges.columns:
        return {"positive": 0, "negative": 0, "unknown": 0, "signed_total": 0}

    positive = 0
    negative = 0
    unknown = 0
    for rel in edges["relation_type"].fillna(""):
        pol = _relation_polarity(str(rel))
        if pol == 1:
            positive += 1
        elif pol == -1:
            negative += 1
        else:
            unknown += 1

    return {
        "positive": positive,
        "negative": negative,
        "unknown": unknown,
        "signed_total": positive + negative,
    }


def _detect_chemical_in_reaction_io(reaction_io: pd.DataFrame) -> bool:
    if reaction_io.empty:
        return False

    type_cols = [c for c in reaction_io.columns if "type" in c.lower()]
    label_cols = [c for c in reaction_io.columns if "label" in c.lower()]

    for col in type_cols:
        if reaction_io[col].fillna("").astype(str).str.contains(
            "SimpleEntity|Chemical|compound", case=False, regex=True
        ).any():
            return True

    chemical_keywords = ["O2", "2OG", "CO2", "SUCCA", "ATP", "ADP", "H2O", "oxygen", "succinate"]
    for col in label_cols:
        text = reaction_io[col].fillna("").astype(str)
        for kw in chemical_keywords:
            if text.str.contains(kw, case=False, regex=False).any():
                return True

    return False


def _target_subgraph_edges(edges: pd.DataFrame, target_node_ids: list[str]) -> pd.DataFrame:
    if edges.empty or not {"source_node", "target_node"}.issubset(edges.columns):
        return pd.DataFrame()
    ids = set(target_node_ids)
    return edges[
        edges["source_node"].astype(str).isin(ids)
        | edges["target_node"].astype(str).isin(ids)
    ]


def evaluate_model_evidence(
    tables: dict[str, Any],
    *,
    target_keyword: str | None = None,
    sbml_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate available mathematical model evidence from loaded result tables."""
    nodes: pd.DataFrame = tables.get("nodes", pd.DataFrame())
    edges: pd.DataFrame = tables.get("edges", pd.DataFrame())
    reaction_io: pd.DataFrame = tables.get("reaction_io", pd.DataFrame())
    event_order: pd.DataFrame = tables.get("event_order_links", pd.DataFrame())
    pathway_events: pd.DataFrame = tables.get("pathway_event_links", pd.DataFrame())
    chemical_nodes: pd.DataFrame = tables.get("chemical_nodes", pd.DataFrame())

    source = ""
    if not nodes.empty and "source" in nodes.columns:
        values = [x for x in nodes["source"].dropna().astype(str).unique().tolist() if x]
        source = ",".join(values[:5])

    target_matches = pd.DataFrame()
    target_node_ids: list[str] = []
    if target_keyword:
        target_matches = classify_targets_by_keyword(nodes, target_keyword)
        if "node_id" in target_matches.columns:
            target_node_ids = target_matches["node_id"].astype(str).tolist()

    target_edges = _target_subgraph_edges(edges, target_node_ids) if target_node_ids else pd.DataFrame()
    signed_all = _signed_edge_stats(edges)
    signed_target = _signed_edge_stats(target_edges)

    node_type_counts = _value_counts(nodes, "type")
    schema_counts = _value_counts(nodes, "schemaClass")
    relation_counts = _value_counts(edges, "relation_type")

    has_graph = bool(len(nodes) > 0 and len(edges) > 0)
    has_signed = bool(signed_all["signed_total"] > 0)
    has_target_signed = bool(signed_target["signed_total"] > 0)

    reaction_node_count = 0
    if "type" in nodes.columns:
        reaction_node_count += int(nodes["type"].fillna("").astype(str).str.contains(
            "Reaction|BlackBoxEvent|ReactionLikeEvent", case=False, regex=True
        ).sum())
    if "schemaClass" in nodes.columns:
        reaction_node_count += int(nodes["schemaClass"].fillna("").astype(str).str.contains(
            "Reaction|BlackBoxEvent|ReactionLikeEvent", case=False, regex=True
        ).sum())

    has_reaction_io = bool(len(reaction_io) > 0 or (
        not edges.empty and "relation_type" in edges.columns
        and edges["relation_type"].fillna("").astype(str).isin(["input", "output"]).any()
    ))

    has_event_order = bool(len(event_order) > 0 or (
        not edges.empty and "relation_type" in edges.columns
        and edges["relation_type"].fillna("").astype(str).eq("precedingEvent").any()
    ))

    has_pathway_event = bool(len(pathway_events) > 0 or (
        not edges.empty and "relation_type" in edges.columns
        and edges["relation_type"].fillna("").astype(str).eq("has_event").any()
    ))

    has_chemicals = bool(len(chemical_nodes) > 0 or (
        "type" in nodes.columns and nodes["type"].fillna("").astype(str).str.contains(
            "compound|SimpleEntity|Chemical", case=False, regex=True
        ).any()
    ))

    stoich_like = bool(has_reaction_io and (len(chemical_nodes) > 0 or _detect_chemical_in_reaction_io(reaction_io)))

    sbml_summary = sbml_summary or {}
    sbml_present = bool(sbml_summary.get("sbml_file_count", 0) or sbml_summary.get("reaction_count", 0))
    kinetic_laws = int(sbml_summary.get("kinetic_law_count", 0) or 0)
    parameter_count = int(sbml_summary.get("parameter_count", 0) or 0)
    rules_count = int(sbml_summary.get("rule_count", 0) or 0)

    kinetic_ode_evidence = bool(sbml_present and (kinetic_laws > 0 or rules_count > 0) and parameter_count > 0)
    sbml_candidate = bool(has_reaction_io or reaction_node_count > 0 or sbml_present)

    # Candidate model types are intentionally conservative.
    model_types = []
    if has_graph:
        model_types.append("Directed labeled graph / topology model")
    if has_signed:
        model_types.append("Signed regulatory network candidate")
    if has_target_signed and target_keyword:
        model_types.append(f"Target-specific signed network candidate: {target_keyword}")
    if has_signed:
        model_types.append("Boolean model candidate")
    if has_reaction_io or reaction_node_count > 0:
        model_types.append("Reaction-event graph model")
    if stoich_like:
        model_types.append("Stoichiometric-like reaction table candidate")
    if sbml_present:
        model_types.append("SBML structural model")
    if kinetic_ode_evidence:
        model_types.append("Kinetic/ODE model evidence")

    conclusion = []
    if has_graph:
        conclusion.append("Parsed pathway can be treated as a graph model G=(V,E).")
    if has_signed:
        conclusion.append("Activation/inhibition/expression edges can be converted into a signed regulatory network.")
    if has_reaction_io:
        conclusion.append("Reaction input/output evidence supports a Reactome-style reaction-event model.")
    if stoich_like:
        conclusion.append("Chemical input/output evidence supports a stoichiometric-like table, but mass balance is not guaranteed.")
    if not kinetic_ode_evidence:
        conclusion.append("No executable ODE/kinetic model is confirmed without SBML kinetic laws and parameters.")
    else:
        conclusion.append("SBML kinetic laws/parameters suggest executable kinetic model evidence.")

    return {
        "source": source,
        "target_keyword": target_keyword,
        "target_matches": target_matches.fillna("").to_dict(orient="records") if not target_matches.empty else [],
        "target_node_ids": target_node_ids,
        "node_count": int(len(nodes)),
        "edge_count": int(len(edges)),
        "node_type_counts": node_type_counts,
        "schema_class_counts": schema_counts,
        "relation_type_counts": relation_counts,
        "signed_edge_stats_all": signed_all,
        "signed_edge_stats_target": signed_target,
        "reaction_io_rows": int(len(reaction_io)),
        "reaction_node_count": int(reaction_node_count),
        "event_order_rows": int(len(event_order)),
        "pathway_event_rows": int(len(pathway_events)),
        "chemical_node_count": int(len(chemical_nodes)),
        "evidence_flags": {
            "topology_graph_model": has_graph,
            "signed_regulatory_model_candidate": has_signed,
            "target_signed_regulatory_model_candidate": has_target_signed,
            "boolean_model_candidate": has_signed,
            "reactome_reaction_event_model": bool(has_reaction_io or reaction_node_count > 0),
            "event_order_model": has_event_order,
            "pathway_event_model": has_pathway_event,
            "stoichiometric_like_model_candidate": stoich_like,
            "sbml_candidate": sbml_candidate,
            "sbml_structural_model_present": sbml_present,
            "kinetic_ode_model_evidence": kinetic_ode_evidence,
        },
        "model_types": model_types,
        "sbml_summary": sbml_summary,
        "conclusion": conclusion,
    }


def evaluate_targets_from_result_dir(
    result_dir: str | Path,
    targets: list[str] | None = None,
    sbml_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate model evidence for one result directory and optional target list."""
    tables = load_result_tables(result_dir)
    if not targets:
        return {"overall": evaluate_model_evidence(tables, sbml_summary=sbml_summary)}

    return {
        target: evaluate_model_evidence(tables, target_keyword=target, sbml_summary=sbml_summary)
        for target in targets
    }


def _md_counts(title: str, counts: dict[str, int]) -> list[str]:
    lines = [f"### {title}", "", "| item | count |", "|---|---:|"]
    if not counts:
        lines.append("| none | 0 |")
        return lines
    for k, v in counts.items():
        lines.append(f"| {k} | {v} |")
    return lines


def write_model_evidence_report(evidence: dict[str, Any], output_md: str | Path) -> None:
    """Write a Markdown report from model evidence dictionary."""
    path = Path(output_md)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = ["# Model Evidence Report", ""]
    items = evidence.items() if any(isinstance(v, dict) and "evidence_flags" in v for v in evidence.values()) else [("overall", evidence)]

    for target, ev in items:
        if not isinstance(ev, dict):
            continue
        lines.append(f"## {target}")
        lines.append("")
        lines.append(f"- source: {ev.get('source', '')}")
        lines.append(f"- target_keyword: {ev.get('target_keyword', '')}")
        lines.append(f"- nodes: {ev.get('node_count', 0)}")
        lines.append(f"- edges: {ev.get('edge_count', 0)}")
        lines.append("")
        lines.extend(_md_counts("Model/evidence flags", {k: int(bool(v)) for k, v in ev.get("evidence_flags", {}).items()}))
        lines.append("")
        lines.append("### Model types")
        for mt in ev.get("model_types", []):
            lines.append(f"- {mt}")
        if not ev.get("model_types"):
            lines.append("- None detected")
        lines.append("")
        lines.append("### Conclusion")
        for c in ev.get("conclusion", []):
            lines.append(f"- {c}")
        lines.append("")

        if ev.get("target_matches"):
            lines.append("### Target matches")
            records = ev["target_matches"][:20]
            keys = list(records[0].keys())
            lines.append("| " + " | ".join(keys) + " |")
            lines.append("|" + "|".join(["---"] * len(keys)) + "|")
            for rec in records:
                vals = [str(rec.get(k, "")).replace("\n", " ")[:120] for k in keys]
                lines.append("| " + " | ".join(vals) + " |")
            lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Evaluate mathematical model evidence from S2SignalPathway results.")
    parser.add_argument("result_dir", help="Result directory containing nodes.csv/edges.csv")
    parser.add_argument("--targets", nargs="*", default=None, help="Target keywords, e.g. HIF1A VEGFA")
    parser.add_argument("--sbml-summary-json", default=None, help="Optional SBML summary JSON")
    parser.add_argument("--output-json", default=None)
    parser.add_argument("--output-md", default=None)
    args = parser.parse_args()

    sbml_summary = {}
    if args.sbml_summary_json:
        sbml_summary = _read_json(Path(args.sbml_summary_json))

    evidence = evaluate_targets_from_result_dir(args.result_dir, args.targets, sbml_summary=sbml_summary)

    if args.output_json:
        Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output_json).write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"Saved: {args.output_json}")
    else:
        print(json.dumps(evidence, ensure_ascii=False, indent=2))

    if args.output_md:
        write_model_evidence_report(evidence, args.output_md)
        print(f"Saved: {args.output_md}")


if __name__ == "__main__":
    main()
