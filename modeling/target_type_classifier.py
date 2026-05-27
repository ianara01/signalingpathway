# -*- coding: utf-8 -*-
"""Target type classifier.

This module classifies parsed KEGG/Reactome nodes into target categories such as
Protein, DNA/Gene, RNA, Complex, Chemical, Reaction/Event, and Pathway.

It is intentionally heuristic because KEGG and Reactome expose different
schemas:
- KEGG: type=gene/compound/map/group
- Reactome: schemaClass/type + reference_database + label/compartment
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class TargetClassification:
    node_id: str
    label: str
    source: str
    node_type: str
    schema_class: str
    target_type: str
    confidence: float
    evidence: str
    reference_database: str = ""
    reference_identifier: str = ""
    compartment: str = ""
    species: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _s(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _lower(value: Any) -> str:
    return _s(value).lower()


def _row_get(row: Any, key: str, default: str = "") -> str:
    """Works for pandas Series, dict, and node-attribute dictionaries."""
    try:
        value = row.get(key, default)
    except AttributeError:
        value = default
    return _s(value)


def infer_target_type_from_row(row: Any, node_id: str | None = None) -> TargetClassification:
    """Infer target type from one node row or node-attribute dictionary."""
    node_id = _s(node_id or _row_get(row, "node_id"))
    label = _row_get(row, "label")
    source = _lower(_row_get(row, "source"))
    node_type = _row_get(row, "type")
    schema_class = _row_get(row, "schemaClass") or node_type
    reference_database = _row_get(row, "reference_database")
    reference_identifier = _row_get(row, "reference_identifier")
    compartment = _row_get(row, "compartment")
    species = _row_get(row, "species")

    node_type_l = node_type.lower()
    schema_l = schema_class.lower()
    label_l = label.lower()
    refdb_l = reference_database.lower()

    # Reactome / reference-based classifications
    if "chebi" in refdb_l or schema_l in {"simpleentity", "chemicaldrug", "chemicalcompound"}:
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            "Chemical", 0.95, "ChEBI/reference database or SimpleEntity/Chemical schema",
            reference_database, reference_identifier, compartment, species,
        )

    if "uniprot" in refdb_l:
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            "Protein", 0.95, "UniProt reference database",
            reference_database, reference_identifier, compartment, species,
        )

    if "ensembl" in refdb_l:
        # Reactome represents genes as EntityWithAccessionedSequence too.
        if "gene" in label_l or reference_identifier.upper().startswith("ENSG"):
            target_type = "DNA/Gene"
            evidence = "ENSEMBL gene identifier or gene label"
        elif reference_identifier.upper().startswith(("ENST", "ENSR")) or "rna" in label_l:
            target_type = "RNA"
            evidence = "ENSEMBL transcript/RNA-like identifier or RNA label"
        else:
            target_type = "DNA/Gene"
            evidence = "ENSEMBL reference database"
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            target_type, 0.9, evidence,
            reference_database, reference_identifier, compartment, species,
        )

    if schema_l == "complex":
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            "Complex", 0.95, "Reactome Complex schema",
            reference_database, reference_identifier, compartment, species,
        )

    if schema_l in {"definedset", "candidateset", "openset"}:
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            "Set/Family", 0.85, "Reactome set schema",
            reference_database, reference_identifier, compartment, species,
        )

    if schema_l in {
        "reaction", "reactionlikeevent", "blackboxevent",
        "polymerisation", "depolymerisation", "failedreaction",
    }:
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            "Reaction/Event", 0.95, "Reactome reaction/event schema",
            reference_database, reference_identifier, compartment, species,
        )

    if schema_l in {"pathway", "toplevelpathway"} or node_type_l == "map":
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            "Pathway", 0.95, "Pathway/map node type",
            reference_database, reference_identifier, compartment, species,
        )

    # KEGG classifications
    if source == "kegg" or node_id.startswith("KEGG:"):
        if node_type_l == "gene":
            # KEGG gene entries usually represent genes/protein products, not explicit RNA/protein.
            return TargetClassification(
                node_id, label, source, node_type, schema_class,
                "Gene/ProteinProduct", 0.85, "KEGG gene entry",
                reference_database, reference_identifier, compartment, species,
            )
        if node_type_l == "compound":
            return TargetClassification(
                node_id, label, source, node_type, schema_class,
                "Chemical", 0.9, "KEGG compound entry",
                reference_database, reference_identifier, compartment, species,
            )
        if node_type_l == "group":
            return TargetClassification(
                node_id, label, source, node_type, schema_class,
                "Group/ComplexLike", 0.75, "KEGG group entry",
                reference_database, reference_identifier, compartment, species,
            )

    # label-based fallback
    if "rna" in label_l or "mirna" in label_l or "lncrna" in label_l:
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            "RNA", 0.6, "RNA keyword in label",
            reference_database, reference_identifier, compartment, species,
        )

    if "gene" in label_l:
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            "DNA/Gene", 0.6, "gene keyword in label",
            reference_database, reference_identifier, compartment, species,
        )

    if "protein" in label_l:
        return TargetClassification(
            node_id, label, source, node_type, schema_class,
            "Protein", 0.55, "protein keyword in label",
            reference_database, reference_identifier, compartment, species,
        )

    return TargetClassification(
        node_id, label, source, node_type, schema_class,
        "Unknown", 0.1, "No reliable schema/reference evidence",
        reference_database, reference_identifier, compartment, species,
    )


def classify_node(node_id: str, attrs: dict[str, Any]) -> dict[str, Any]:
    """Classify one NetworkX node using its attributes."""
    data = dict(attrs)
    data["node_id"] = node_id
    return infer_target_type_from_row(data).to_dict()


def classify_nodes_dataframe(nodes: pd.DataFrame) -> pd.DataFrame:
    """Add target classification for every row of nodes.csv."""
    if nodes.empty:
        return pd.DataFrame()

    records = []
    for _, row in nodes.iterrows():
        records.append(infer_target_type_from_row(row).to_dict())

    return pd.DataFrame(records)


def classify_targets_by_keyword(
    nodes: pd.DataFrame,
    keyword: str,
    *,
    search_columns: list[str] | None = None,
    max_rows: int = 100,
) -> pd.DataFrame:
    """Find nodes matching keyword and classify them."""
    if nodes.empty:
        return pd.DataFrame()

    if search_columns is None:
        search_columns = [
            "node_id", "label", "raw_id", "reference_identifier",
            "reference_display", "compound_name", "compound_id",
        ]

    available = [c for c in search_columns if c in nodes.columns]
    if not available:
        return pd.DataFrame()

    mask = pd.Series(False, index=nodes.index)
    for col in available:
        mask = mask | nodes[col].fillna("").astype(str).str.contains(
            keyword, case=False, regex=False
        )

    matched = nodes[mask].head(max_rows)
    return classify_nodes_dataframe(matched)


def load_nodes_csv(path: str | Path) -> pd.DataFrame:
    """Load nodes.csv with encoding fallback."""
    p = Path(path)
    try:
        return pd.read_csv(p, encoding="utf-8-sig")
    except UnicodeDecodeError:
        return pd.read_csv(p)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Classify KEGG/Reactome target nodes.")
    parser.add_argument("nodes_csv", help="Path to nodes.csv")
    parser.add_argument("--keyword", default=None, help="Optional target keyword, e.g. HIF1A")
    parser.add_argument("--output", default=None, help="Output CSV path")
    args = parser.parse_args()

    nodes = load_nodes_csv(args.nodes_csv)
    if args.keyword:
        out = classify_targets_by_keyword(nodes, args.keyword)
    else:
        out = classify_nodes_dataframe(nodes)

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.output, index=False, encoding="utf-8-sig")
        print(f"Saved: {args.output}")
    else:
        print(out.to_string(index=False))


if __name__ == "__main__":
    main()
