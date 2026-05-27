# -*- coding: utf-8 -*-
"""Modeling utilities for S2SignalPathway.

This package does not replace KEGG/Reactome parsers. It inspects parsed
nodes/edges/reaction tables and reports what kind of mathematical model can be
constructed from the available evidence.

Main model levels:
- graph/topology model
- signed regulatory network
- Boolean rule candidate
- Reactome reaction-event graph
- stoichiometric-like reaction table
- SBML/kinetic model evidence
"""

from .target_type_classifier import (
    classify_node,
    classify_nodes_dataframe,
    classify_targets_by_keyword,
    infer_target_type_from_row,
)

from .model_evidence_checker import (
    load_result_tables,
    evaluate_model_evidence,
    evaluate_targets_from_result_dir,
    write_model_evidence_report,
)

from .sbml_checker import (
    summarize_sbml,
    find_sbml_files,
    summarize_sbml_directory,
)

from .boolean_model_builder import (
    build_signed_regulatory_edges,
    build_boolean_rules,
    write_bnet_rules,
    write_boolean_rule_report,
)

__all__ = [
    "classify_node",
    "classify_nodes_dataframe",
    "classify_targets_by_keyword",
    "infer_target_type_from_row",
    "load_result_tables",
    "evaluate_model_evidence",
    "evaluate_targets_from_result_dir",
    "write_model_evidence_report",
    "summarize_sbml",
    "find_sbml_files",
    "summarize_sbml_directory",
    "build_signed_regulatory_edges",
    "build_boolean_rules",
    "write_bnet_rules",
    "write_boolean_rule_report",
]
