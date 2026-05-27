# -*- coding: utf-8 -*-
"""SBML checker.

This module performs a lightweight SBML structural/kinetic-evidence check using
Python's standard XML parser. It does not require libSBML.

It checks:
- species
- reactions
- parameters
- compartments
- kineticLaw elements
- rules/events
- initialAmount/initialConcentration presence

Interpretation:
- SBML structural model exists if species and reactions exist.
- kinetic/ODE evidence is stronger if kineticLaw + parameters/rules are present.
- Full simulation validity still requires a proper SBML validator/libSBML.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import xml.etree.ElementTree as ET


def _strip_ns(tag: str) -> str:
    return tag.split("}", 1)[-1] if "}" in tag else tag


def _find_all_by_local_name(root: ET.Element, local_name: str) -> list[ET.Element]:
    return [el for el in root.iter() if _strip_ns(el.tag) == local_name]


def _attr(el: ET.Element, *names: str) -> str:
    for name in names:
        if name in el.attrib:
            return str(el.attrib[name])
    return ""


def summarize_sbml(sbml_path: str | Path, *, max_items: int = 20) -> dict[str, Any]:
    """Summarize one SBML file."""
    path = Path(sbml_path)
    if not path.exists():
        raise FileNotFoundError(f"SBML file not found: {path}")

    tree = ET.parse(path)
    root = tree.getroot()

    model_nodes = _find_all_by_local_name(root, "model")
    model = model_nodes[0] if model_nodes else None

    species = _find_all_by_local_name(root, "species")
    reactions = _find_all_by_local_name(root, "reaction")
    parameters = _find_all_by_local_name(root, "parameter")
    compartments = _find_all_by_local_name(root, "compartment")
    kinetic_laws = _find_all_by_local_name(root, "kineticLaw")
    rules = (
        _find_all_by_local_name(root, "assignmentRule")
        + _find_all_by_local_name(root, "rateRule")
        + _find_all_by_local_name(root, "algebraicRule")
    )
    events = _find_all_by_local_name(root, "event")

    initial_value_count = 0
    species_samples = []
    for sp in species:
        if "initialAmount" in sp.attrib or "initialConcentration" in sp.attrib:
            initial_value_count += 1
        if len(species_samples) < max_items:
            species_samples.append({
                "id": _attr(sp, "id"),
                "name": _attr(sp, "name"),
                "compartment": _attr(sp, "compartment"),
                "initialAmount": _attr(sp, "initialAmount"),
                "initialConcentration": _attr(sp, "initialConcentration"),
                "boundaryCondition": _attr(sp, "boundaryCondition"),
                "constant": _attr(sp, "constant"),
            })

    reaction_samples = []
    for rxn in reactions[:max_items]:
        reaction_samples.append({
            "id": _attr(rxn, "id"),
            "name": _attr(rxn, "name"),
            "reversible": _attr(rxn, "reversible"),
            "has_kineticLaw": any(_strip_ns(child.tag) == "kineticLaw" for child in list(rxn)),
        })

    parameter_samples = []
    for par in parameters[:max_items]:
        parameter_samples.append({
            "id": _attr(par, "id"),
            "name": _attr(par, "name"),
            "value": _attr(par, "value"),
            "units": _attr(par, "units"),
            "constant": _attr(par, "constant"),
        })

    species_count = len(species)
    reaction_count = len(reactions)
    parameter_count = len(parameters)
    kinetic_law_count = len(kinetic_laws)
    rule_count = len(rules)
    event_count = len(events)

    structural_model = bool(species_count > 0 and reaction_count > 0)
    kinetic_evidence = bool(kinetic_law_count > 0 or rule_count > 0)
    parameterized = bool(parameter_count > 0)
    executable_candidate = bool(structural_model and kinetic_evidence and (parameterized or initial_value_count > 0))

    limitations = []
    if not structural_model:
        limitations.append("No species/reaction structure detected.")
    if kinetic_law_count == 0:
        limitations.append("No kineticLaw elements detected.")
    if parameter_count == 0:
        limitations.append("No SBML parameter elements detected.")
    if initial_value_count == 0:
        limitations.append("No species initialAmount/initialConcentration values detected.")
    limitations.append("This is a lightweight XML check; use libSBML for validation and simulation readiness.")

    return {
        "file": str(path),
        "model_id": _attr(model, "id") if model is not None else "",
        "model_name": _attr(model, "name") if model is not None else "",
        "species_count": species_count,
        "reaction_count": reaction_count,
        "parameter_count": parameter_count,
        "compartment_count": len(compartments),
        "kinetic_law_count": kinetic_law_count,
        "rule_count": rule_count,
        "event_count": event_count,
        "species_initial_value_count": initial_value_count,
        "has_structural_model": structural_model,
        "has_kinetic_law_or_rules": kinetic_evidence,
        "has_parameters": parameterized,
        "kinetic_ode_model_candidate": executable_candidate,
        "species_samples": species_samples,
        "reaction_samples": reaction_samples,
        "parameter_samples": parameter_samples,
        "limitations": limitations,
    }


def find_sbml_files(root_dir: str | Path) -> list[Path]:
    """Find .xml/.sbml files below root_dir."""
    root = Path(root_dir)
    if not root.exists():
        return []

    candidates: list[Path] = []
    for pattern in ("*.sbml", "*.xml", "*.SBML", "*.XML"):
        candidates.extend(root.rglob(pattern))

    # Keep files that appear to contain SBML.
    sbml_files: list[Path] = []
    for path in sorted(set(candidates)):
        try:
            head = path.read_text(encoding="utf-8", errors="ignore")[:1000].lower()
            if "<sbml" in head:
                sbml_files.append(path)
        except Exception:
            continue

    return sbml_files


def summarize_sbml_directory(root_dir: str | Path) -> dict[str, Any]:
    """Summarize all SBML files below a directory."""
    files = find_sbml_files(root_dir)
    summaries = []
    for f in files:
        try:
            summaries.append(summarize_sbml(f))
        except Exception as exc:
            summaries.append({"file": str(f), "error": str(exc)})

    return {
        "root_dir": str(root_dir),
        "sbml_file_count": len(files),
        "files": [str(f) for f in files],
        "summaries": summaries,
        "total_species_count": int(sum(x.get("species_count", 0) for x in summaries if isinstance(x, dict))),
        "total_reaction_count": int(sum(x.get("reaction_count", 0) for x in summaries if isinstance(x, dict))),
        "total_parameter_count": int(sum(x.get("parameter_count", 0) for x in summaries if isinstance(x, dict))),
        "total_kinetic_law_count": int(sum(x.get("kinetic_law_count", 0) for x in summaries if isinstance(x, dict))),
        "any_kinetic_ode_model_candidate": any(bool(x.get("kinetic_ode_model_candidate")) for x in summaries if isinstance(x, dict)),
    }


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Lightweight SBML model checker.")
    parser.add_argument("path", help="SBML file or directory")
    parser.add_argument("--output-json", default=None)
    args = parser.parse_args()

    p = Path(args.path)
    if p.is_dir():
        result = summarize_sbml_directory(p)
    else:
        result = summarize_sbml(p)

    if args.output_json:
        Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output_json).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"Saved: {args.output_json}")
    else:
        print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
