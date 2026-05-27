# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:11:17 2026

@author: user
"""

# workflows/target_registry.py
# workflows/target_registry.py

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


TARGET_REGISTRY: dict[str, dict[str, Any]] = {
    "hif1": {
        "aliases": [
            "hif-1",
            "hif1",
            "hif1a",
            "hypoxia",
            "hypoxia inducible factor",
            "hypoxia-inducible factor",
        ],
        "display_name": "HIF-1 / Hypoxia-inducible Factor",
        "kegg": {
            "pathway_id": "hsa04066",
            "pathway_name": "HIF-1 signaling pathway",
            "keywords": [
                "HIF1A", "HIF-1", "VEGF", "VEGFA", "MTOR",
                "ARNT", "VHL", "EGLN", "PHD", "PDK1", "SLC2A1",
            ],
            "representative_paths": [
                {
                    "name": "HIF1A_to_VEGFA",
                    "start": "KEGG:3091",
                    "end": "KEGG:7422",
                }
            ],
        },
        "reactome": {
            "steps": [
                {
                    "name": "01_pathway_cellular_response_to_hypoxia",
                    "reactome_id": "R-HSA-1234174",
                    "description": "Cellular response to hypoxia",
                    "build_reactome_flow": True,
                    "recursive_events": True,
                    "recursive_max_depth": 3,
                },
                {
                    "name": "02_entity_hif1a_nucleoplasm",
                    "reactome_id": "R-HSA-1234135",
                    "description": "HIF1A [nucleoplasm]",
                    "build_reactome_flow": True,
                    "recursive_events": False,
                    "recursive_max_depth": 0,
                },
                {
                    "name": "03_reaction_hif_alpha_binds_arnt",
                    "reactome_id": "R-HSA-1234171",
                    "description": "HIF-alpha binds ARNT forming HIF-alpha:ARNT",
                    "build_reactome_flow": True,
                    "recursive_events": False,
                    "recursive_max_depth": 0,
                },
            ],
            "keywords": [
                "HIF1A", "HIF-1", "HIF", "ARNT", "VEGF", "VEGFA",
                "EPO", "CA9", "O2", "2OG", "CO2", "SUCCA", "VHL", "PHD",
            ],
            "representative_paths": [
                {
                    "name": "HIF1A_to_HIF_alpha_binds_ARNT",
                    "start": "REACTOME:R-HSA-1234135",
                    "end": "REACTOME:R-HSA-1234171",
                }
            ],
        },
    },
}


CUSTOM_REGISTRY_PATH = Path("workflows") / "target_registry_custom.json"


def normalize_target_text(text: str) -> str:
    return (
        str(text)
        .strip()
        .lower()
        .replace("_", "-")
        .replace(" ", "")
    )


def load_custom_registry() -> dict[str, dict[str, Any]]:
    if not CUSTOM_REGISTRY_PATH.exists():
        return {}

    with CUSTOM_REGISTRY_PATH.open("r", encoding="utf-8") as f:
        return json.load(f)


def merged_registry() -> dict[str, dict[str, Any]]:
    registry = dict(TARGET_REGISTRY)
    registry.update(load_custom_registry())
    return registry


def recognize_target_key(target_text: str) -> str | None:
    normalized = normalize_target_text(target_text)

    for target_key, config in merged_registry().items():
        aliases = config.get("aliases", [])

        for alias in aliases:
            alias_norm = normalize_target_text(alias)

            if normalized == alias_norm:
                return target_key

        for alias in aliases:
            alias_norm = normalize_target_text(alias)

            if alias_norm and alias_norm in normalized:
                return target_key

    return None


def get_target_config(target_key: str) -> dict[str, Any]:
    registry = merged_registry()
    return registry[target_key]


def save_custom_target(target_key: str, config: dict[str, Any]) -> None:
    custom = load_custom_registry()
    custom[target_key] = config

    CUSTOM_REGISTRY_PATH.parent.mkdir(parents=True, exist_ok=True)

    with CUSTOM_REGISTRY_PATH.open("w", encoding="utf-8") as f:
        json.dump(custom, f, ensure_ascii=False, indent=2)
