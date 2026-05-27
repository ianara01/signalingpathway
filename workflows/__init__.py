# -*- coding: utf-8 -*-
"""
Created on Mon May 25 11:33:23 2026

@author: user
"""

from .target_registry import (
    TARGET_REGISTRY,
    recognize_target_key,
    get_target_config,
    save_custom_target,
)

from .generic_kegg_target_workflow import run_generic_kegg_target_workflow
from .generic_reactome_target_workflow import run_generic_reactome_target_workflow

__all__ = [
    "TARGET_REGISTRY",
    "recognize_target_key",
    "get_target_config",
    "save_custom_target",
    "run_generic_kegg_target_workflow",
    "run_generic_reactome_target_workflow",
]
