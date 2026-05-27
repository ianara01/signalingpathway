# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:11:17 2026

@author: user
"""

# workflows/target_registry.py


TARGET_REGISTRY = {
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
                "HIF1A",
                "HIF-1",
                "VEGF",
                "VEGFA",
                "MTOR",
                "ARNT",
                "VHL",
                "EGLN",
                "PHD",
                "PDK1",
                "SLC2A1",
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
                "HIF1A",
                "HIF-1",
                "HIF",
                "ARNT",
                "VEGF",
                "VEGFA",
                "EPO",
                "CA9",
                "O2",
                "2OG",
                "CO2",
                "SUCCA",
                "VHL",
                "PHD",
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

    "p53": {
        "aliases": [
            "p53",
            "tp53",
            "tumor protein p53",
            "p53 signaling",
        ],
        "display_name": "p53 / TP53 signaling",
        "kegg": {
            "pathway_id": "hsa04115",
            "pathway_name": "p53 signaling pathway",
            "keywords": [
                "TP53",
                "p53",
                "MDM2",
                "CDKN1A",
                "BAX",
                "BBC3",
                "GADD45",
                "ATM",
                "ATR",
                "CHEK1",
                "CHEK2",
                "CASP",
            ],
            "representative_paths": [
                {
                    "name": "TP53_to_CDKN1A",
                    "start": "KEGG:7157",
                    "end": "KEGG:1026",
                },
                {
                    "name": "TP53_to_BAX",
                    "start": "KEGG:7157",
                    "end": "KEGG:581",
                },
            ],
        },
        "reactome": {
            "steps": [],
            "keywords": [
                "TP53",
                "p53",
                "MDM2",
                "CDKN1A",
                "BAX",
                "DNA damage",
                "apoptosis",
                "cell cycle",
            ],
            "representative_paths": [],
        },
    },
}


def normalize_target_text(text: str) -> str:
    return (
        str(text)
        .strip()
        .lower()
        .replace("_", "-")
        .replace(" ", "")
    )

def recognize_target_key(target_text: str) -> str | None:
    normalized = normalize_target_text(target_text)

    for target_key, config in TARGET_REGISTRY.items():
        aliases = config.get("aliases", [])

        for alias in aliases:
            alias_norm = normalize_target_text(alias)

            if normalized == alias_norm:
                return target_key

        # 부분 매칭
        for alias in aliases:
            alias_norm = normalize_target_text(alias)

            if alias_norm and alias_norm in normalized:
                return target_key

    return None

def get_target_config(target_key: str) -> dict:
    return TARGET_REGISTRY[target_key]

