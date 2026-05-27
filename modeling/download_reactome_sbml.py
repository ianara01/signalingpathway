# -*- coding: utf-8 -*-
"""
Created on Wed May 27 14:20:53 2026

@author: user


download_reactome_sbml.py

Reactome pathway/reaction event ID를 SBML로 다운로드하여 data/sbml/reactome/{target_key}/에 저장한다.

사용 예:
  python download_reactome_sbml.py ^
    --target-key hif1 ^
    --ids R-HSA-1234174 R-HSA-1234171 ^
    --output-root data/sbml/reactome
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import requests


REACTOME_BASE_URL = "https://reactome.org/ContentService"


def safe_filename(reactome_id: str) -> str:
    return (
        str(reactome_id)
        .strip()
        .replace(":", "_")
        .replace("/", "_")
        .replace("\\", "_")
    )


def query_reactome_id(
    reactome_id: str,
    *,
    timeout: tuple[int, int] = (10, 120),
) -> dict[str, Any]:
    """
    Reactome ID의 schemaClass를 확인한다.
    """
    url = f"{REACTOME_BASE_URL}/data/query/{reactome_id}"

    response = requests.get(
        url,
        timeout=timeout,
        headers={
            "Accept": "application/json",
            "User-Agent": "S2SignalPathway/1.0",
        },
    )
    response.raise_for_status()
    return response.json()


def is_exportable_event(schema_class: str) -> bool:
    """
    SBML export 대상이 될 수 있는 Reactome event 계열인지 확인한다.
    """
    return schema_class in {
        "Pathway",
        "TopLevelPathway",
        "Reaction",
        "ReactionLikeEvent",
        "BlackBoxEvent",
        "Polymerisation",
        "Depolymerisation",
    }


def download_reactome_event_sbml(
    reactome_id: str,
    output_dir: Path,
    *,
    timeout: tuple[int, int] = (10, 180),
    overwrite: bool = False,
    sleep_sec: float = 0.5,
) -> dict[str, Any]:
    """
    Reactome event/pathway/reaction ID를 SBML로 다운로드한다.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    clean_id = reactome_id.strip()
    output_path = output_dir / f"{safe_filename(clean_id)}.sbml"

    if output_path.exists() and not overwrite:
        return {
            "reactome_id": clean_id,
            "status": "skipped_exists",
            "output_path": str(output_path),
        }

    # 먼저 schemaClass 확인
    try:
        meta = query_reactome_id(clean_id, timeout=timeout)
        schema_class = meta.get("schemaClass", "")
        display_name = meta.get("displayName", "")
    except Exception as exc:
        return {
            "reactome_id": clean_id,
            "status": "query_failed",
            "error": str(exc),
        }

    if not is_exportable_event(schema_class):
        return {
            "reactome_id": clean_id,
            "schemaClass": schema_class,
            "displayName": display_name,
            "status": "not_exportable_entity",
            "message": (
                "이 ID는 Pathway/Reaction/Event 계열이 아니므로 SBML export를 건너뜁니다. "
                "Entity target이면 관련 pathway/reaction event ID를 먼저 선택하세요."
            ),
        }

    url = f"{REACTOME_BASE_URL}/exporter/event/{clean_id}.sbml"

    try:
        response = requests.get(
            url,
            timeout=timeout,
            headers={
                "Accept": "application/xml, text/xml, application/sbml+xml, */*",
                "User-Agent": "S2SignalPathway/1.0",
            },
        )
        response.raise_for_status()

        text = response.text

        # SBML 여부 1차 확인
        if "<sbml" not in text[:2000].lower():
            return {
                "reactome_id": clean_id,
                "schemaClass": schema_class,
                "displayName": display_name,
                "status": "downloaded_but_not_sbml",
                "url": url,
                "preview": text[:500],
            }

        output_path.write_text(text, encoding="utf-8")

        time.sleep(sleep_sec)

        return {
            "reactome_id": clean_id,
            "schemaClass": schema_class,
            "displayName": display_name,
            "status": "success",
            "url": url,
            "output_path": str(output_path),
            "size_bytes": output_path.stat().st_size,
        }

    except Exception as exc:
        return {
            "reactome_id": clean_id,
            "schemaClass": schema_class,
            "displayName": display_name,
            "status": "download_failed",
            "url": url,
            "error": str(exc),
        }


def download_many_reactome_sbml(
    reactome_ids: list[str],
    target_key: str,
    *,
    output_root: str | Path = "data/sbml/reactome",
    overwrite: bool = False,
) -> dict[str, Any]:
    output_dir = Path(output_root) / target_key
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []

    for reactome_id in reactome_ids:
        print(f"Downloading SBML: {reactome_id}")
        result = download_reactome_event_sbml(
            reactome_id,
            output_dir,
            overwrite=overwrite,
        )
        print(f"  -> {result.get('status')}")
        results.append(result)

    summary = {
        "target_key": target_key,
        "output_dir": str(output_dir),
        "requested_ids": reactome_ids,
        "success_count": sum(1 for r in results if r.get("status") == "success"),
        "skipped_count": sum(1 for r in results if r.get("status", "").startswith("skipped")),
        "failed_count": sum(1 for r in results if "failed" in r.get("status", "")),
        "not_exportable_count": sum(1 for r in results if r.get("status") == "not_exportable_entity"),
        "results": results,
    }

    summary_path = output_dir / "download_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(f"\nDownload summary saved: {summary_path}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-key", required=True, help="예: hif1, p53, keratin")
    parser.add_argument(
        "--ids",
        nargs="+",
        required=True,
        help="Reactome event/pathway/reaction IDs, e.g. R-HSA-1234174 R-HSA-1234171",
    )
    parser.add_argument(
        "--output-root",
        default="data/sbml/reactome",
        help="SBML 저장 루트 디렉토리",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    download_many_reactome_sbml(
        reactome_ids=args.ids,
        target_key=args.target_key,
        output_root=args.output_root,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()