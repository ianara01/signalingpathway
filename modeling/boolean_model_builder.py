# -*- coding: utf-8 -*-
"""Boolean model builder.

This module converts KEGG/Reactome labeled edges into a simple signed network
and generates Boolean-rule candidates.

Important limitations:
- Rules are topology-derived candidates, not validated kinetic models.
- via_group edges are inferred parser expansions and can be excluded by default.
- "expression" is treated as positive regulatory evidence by default.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import json
import re
import pandas as pd


@dataclass
class SignedEdge:
    source_node: str
    target_node: str
    relation_type: str
    sign: int
    confidence: float
    evidence: str
    via_group: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _read_csv(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    try:
        return pd.read_csv(p, encoding="utf-8-sig")
    except UnicodeDecodeError:
        return pd.read_csv(p)


def infer_edge_sign(
    relation_type: str,
    *,
    expression_positive: bool = True,
) -> tuple[int | None, float, str]:
    """Infer sign from relation type."""
    r = str(relation_type).lower()

    if any(x in r for x in ["inhibition", "repression", "negative", "ubiquitination"]):
        return -1, 0.85, "negative relation keyword"

    if "activation" in r or "positive" in r:
        return 1, 0.85, "positive relation keyword"

    if "expression" in r:
        if expression_positive:
            return 1, 0.65, "expression treated as positive regulatory evidence"
        return None, 0.3, "expression relation left unsigned"

    if "phosphorylation" in r:
        return 1, 0.55, "phosphorylation treated as weak positive evidence"

    return None, 0.1, "no sign keyword"


def build_signed_regulatory_edges(
    edges: pd.DataFrame,
    *,
    include_via_group: bool = False,
    expression_positive: bool = True,
) -> pd.DataFrame:
    """Build signed regulatory edge table from edges.csv."""
    if edges.empty or not {"source_node", "target_node"}.issubset(edges.columns):
        return pd.DataFrame()

    records = []
    for _, row in edges.iterrows():
        rel = str(row.get("relation_type", ""))
        via_group = "via_group" in rel.lower() or str(row.get("expanded_from_group", "")).lower() == "true"

        if via_group and not include_via_group:
            continue

        sign, confidence, evidence = infer_edge_sign(rel, expression_positive=expression_positive)
        if sign is None:
            continue

        records.append(SignedEdge(
            source_node=str(row.get("source_node", "")),
            target_node=str(row.get("target_node", "")),
            relation_type=rel,
            sign=sign,
            confidence=confidence,
            evidence=evidence,
            via_group=via_group,
        ).to_dict())

    return pd.DataFrame(records)


def _node_symbol(node_id: str, label: str | None = None) -> str:
    """Create a Boolean-safe symbol from label or node_id."""
    base = (label or "").strip() or str(node_id)
    # Prefer concise gene-like label first token.
    base = base.split("[")[0].strip()
    base = base.split(",")[0].strip()
    base = base.replace("TITLE:", "")
    if not base:
        base = str(node_id)

    base = re.sub(r"[^0-9A-Za-z_]+", "_", base)
    base = re.sub(r"_+", "_", base).strip("_")
    if not base:
        base = re.sub(r"[^0-9A-Za-z_]+", "_", str(node_id))
    if base and base[0].isdigit():
        base = "N_" + base
    return base


def _build_label_map(nodes: pd.DataFrame) -> dict[str, str]:
    if nodes.empty or "node_id" not in nodes.columns:
        return {}
    label_map = {}
    for _, row in nodes.iterrows():
        node_id = str(row.get("node_id", ""))
        label = str(row.get("label", "")) if "label" in nodes.columns else ""
        label_map[node_id] = _node_symbol(node_id, label)
    return label_map


def build_boolean_rules(
    nodes: pd.DataFrame,
    signed_edges: pd.DataFrame,
    *,
    only_targets_with_inputs: bool = True,
) -> dict[str, Any]:
    """Generate simple Boolean rules from signed edges.

    Rule structure:
      target = (activator1 or activator2 ...) and not (inhibitor1 or inhibitor2 ...)

    Nodes with no signed inputs are omitted by default.
    """
    label_map = _build_label_map(nodes)
    if signed_edges.empty:
        return {"rules": {}, "records": [], "limitations": ["No signed edges available."]}

    rules: dict[str, str] = {}
    records: list[dict[str, Any]] = []

    for target_node, group in signed_edges.groupby("target_node"):
        positives = []
        negatives = []
        for _, row in group.iterrows():
            source = str(row["source_node"])
            sign = int(row["sign"])
            symbol = label_map.get(source, _node_symbol(source))
            if sign > 0:
                positives.append(symbol)
            elif sign < 0:
                negatives.append(symbol)

        target_symbol = label_map.get(str(target_node), _node_symbol(str(target_node)))

        pos_expr = " or ".join(sorted(set(positives)))
        neg_expr = " or ".join(sorted(set(negatives)))

        if pos_expr and neg_expr:
            expr = f"({pos_expr}) and not ({neg_expr})"
        elif pos_expr:
            expr = f"({pos_expr})"
        elif neg_expr:
            expr = f"not ({neg_expr})"
        else:
            if only_targets_with_inputs:
                continue
            expr = target_symbol

        rules[target_symbol] = expr
        records.append({
            "target_node": target_node,
            "target_symbol": target_symbol,
            "rule": expr,
            "positive_inputs": ";".join(sorted(set(positives))),
            "negative_inputs": ";".join(sorted(set(negatives))),
            "input_edge_count": int(len(group)),
        })

    limitations = [
        "Rules are topology-derived candidates, not validated dynamic models.",
        "Logical operators are inferred from relation labels only.",
        "No time delay, rate, threshold, or parameter values are included.",
    ]

    return {"rules": rules, "records": records, "limitations": limitations}


def write_bnet_rules(rules: dict[str, str], output_path: str | Path) -> None:
    """Write rules in simple .bnet-like format: target, rule."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["targets, factors"]
    for target, expr in rules.items():
        lines.append(f"{target}, {expr}")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_boolean_rule_report(result: dict[str, Any], output_md: str | Path) -> None:
    path = Path(output_md)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = ["# Boolean Rule Candidate Report", ""]
    lines.append("## Rules")
    lines.append("")
    if not result.get("records"):
        lines.append("- No rules generated.")
    else:
        keys = ["target_symbol", "rule", "positive_inputs", "negative_inputs", "input_edge_count"]
        lines.append("| " + " | ".join(keys) + " |")
        lines.append("|" + "|".join(["---"] * len(keys)) + "|")
        for rec in result["records"]:
            lines.append("| " + " | ".join(str(rec.get(k, "")) for k in keys) + " |")
    lines.append("")
    lines.append("## Limitations")
    for lim in result.get("limitations", []):
        lines.append(f"- {lim}")

    path.write_text("\n".join(lines), encoding="utf-8")


def build_boolean_model_from_csv(
    nodes_csv: str | Path,
    edges_csv: str | Path,
    *,
    include_via_group: bool = False,
    expression_positive: bool = True,
) -> dict[str, Any]:
    nodes = _read_csv(nodes_csv)
    edges = _read_csv(edges_csv)

    signed = build_signed_regulatory_edges(
        edges,
        include_via_group=include_via_group,
        expression_positive=expression_positive,
    )
    result = build_boolean_rules(nodes, signed)
    result["signed_edges"] = signed.fillna("").to_dict(orient="records")
    result["include_via_group"] = include_via_group
    result["expression_positive"] = expression_positive
    return result


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Build Boolean-rule candidates from S2SignalPathway nodes/edges CSV.")
    parser.add_argument("--nodes", required=True, help="nodes.csv")
    parser.add_argument("--edges", required=True, help="edges.csv")
    parser.add_argument("--include-via-group", action="store_true", help="Include KEGG via_group inferred edges")
    parser.add_argument("--expression-unsigned", action="store_true", help="Do not treat expression as positive")
    parser.add_argument("--output-json", default=None)
    parser.add_argument("--output-bnet", default=None)
    parser.add_argument("--output-md", default=None)
    args = parser.parse_args()

    result = build_boolean_model_from_csv(
        args.nodes,
        args.edges,
        include_via_group=args.include_via_group,
        expression_positive=not args.expression_unsigned,
    )

    if args.output_json:
        Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output_json).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"Saved: {args.output_json}")

    if args.output_bnet:
        write_bnet_rules(result.get("rules", {}), args.output_bnet)
        print(f"Saved: {args.output_bnet}")

    if args.output_md:
        write_boolean_rule_report(result, args.output_md)
        print(f"Saved: {args.output_md}")

    if not any([args.output_json, args.output_bnet, args.output_md]):
        print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
