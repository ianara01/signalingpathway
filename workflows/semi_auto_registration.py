# -*- coding: utf-8 -*-
"""
Created on Mon May 25 17:59:07 2026

@author: user
"""

# workflows/semi_auto_registration.py

from __future__ import annotations

from typing import Any

from SignalPathway.s2pathway import SignalingPathway
from workflows.target_registry import (
    normalize_target_text,
    save_custom_target,
)


def _safe_target_key(target_text: str) -> str:
    key = normalize_target_text(target_text)
    key = key.replace("-", "_")
    return key or "custom_target"


def _search_kegg_candidates(target_text: str) -> list[dict[str, Any]]:
    pathway = SignalingPathway(source="kegg")

    try:
        results = pathway.search_pathway_id(target_text)
    except Exception:
        return []

    return results or []


def _select_kegg_candidate(candidates: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not candidates:
        print("\nKEGG 후보가 없습니다.")
        return None

    print("\n--- KEGG 후보 ---")
    for index, item in enumerate(candidates):
        print(f"[{index}] {item.get('id')}: {item.get('name')}")

    text = input("사용할 KEGG 후보 index 또는 pathway ID 입력, skip은 Enter: ").strip()

    if not text:
        return None

    if text.isdigit():
        idx = int(text)
        if 0 <= idx < len(candidates):
            return candidates[idx]
        return None

    for item in candidates:
        if text.lower() == str(item.get("id", "")).lower():
            return item

    if text.lower().startswith(("hsa", "map")):
        return {"id": text, "name": ""}

    return None


def _manual_reactome_steps() -> list[dict[str, Any]]:
    print("\nReactome 후보 자동 검색은 아직 제한적입니다.")
    print("Reactome ID를 알고 있으면 직접 등록할 수 있습니다.")
    print("예: R-HSA-1234174, R-HSA-1234135, R-HSA-1234171")
    print("입력을 끝내려면 빈 줄을 입력하세요.")

    steps = []

    while True:
        reactome_id = input("Reactome ID 입력, 종료는 Enter: ").strip()

        if not reactome_id:
            break

        description = input("설명 입력: ").strip()
        mode = input("유형 [1: pathway, 2: entity, 3: reaction], 기본 1: ").strip()

        if mode == "1":
            recursive_events = input("containedEvents 재귀 확장? (y/n, 기본 n): ").strip().lower() == "y"
            depth = 3
            if recursive_events:
                depth_text = input("recursive depth 기본 3: ").strip()
                if depth_text:
                    depth = int(depth_text)

            steps.append({
                "name": f"{len(steps)+1:02d}_{reactome_id.replace(':', '_')}",
                "reactome_id": reactome_id,
                "description": description,
                "build_reactome_flow": True,
                "recursive_events": recursive_events,
                "recursive_max_depth": depth,
            })

        elif mode == "2":
            steps.append({
                "name": f"{len(steps)+1:02d}_{reactome_id.replace(':', '_')}",
                "reactome_id": reactome_id,
                "description": description,
                "build_reactome_flow": True,
                "recursive_events": False,
                "recursive_max_depth": 0,
            })

        else:
            steps.append({
                "name": f"{len(steps)+1:02d}_{reactome_id.replace(':', '_')}",
                "reactome_id": reactome_id,
                "description": description,
                "build_reactome_flow": True,
                "recursive_events": False,
                "recursive_max_depth": 0,
            })

    return steps


def semi_auto_register_target(target_text: str) -> str | None:
    print(f"\n--- Semi-auto registration: {target_text} ---")

    target_key = _safe_target_key(target_text)
    display_name = input(f"표시 이름 입력, 기본 '{target_text}': ").strip() or target_text

    aliases_text = input("alias들을 쉼표로 입력, 기본 target 포함: ").strip()
    aliases = [target_text]

    if aliases_text:
        aliases.extend([x.strip() for x in aliases_text.split(",") if x.strip()])

    # KEGG 후보 검색
    kegg_candidates = _search_kegg_candidates(target_text)
    kegg_choice = _select_kegg_candidate(kegg_candidates)

    kegg_config = None

    if kegg_choice:
        pathway_id = kegg_choice.get("id", "")
        pathway_name = kegg_choice.get("name", "")

        keywords_text = input(
            "KEGG keyword들을 쉼표로 입력, 기본 target만 사용: "
        ).strip()

        keywords = [target_text]
        if keywords_text:
            keywords.extend([x.strip() for x in keywords_text.split(",") if x.strip()])

        kegg_config = {
            "pathway_id": pathway_id.replace("map", "hsa", 1) if pathway_id.startswith("map") else pathway_id,
            "pathway_name": pathway_name,
            "keywords": keywords,
            "representative_paths": [],
        }

    # Reactome steps 수동/반자동 등록
    add_reactome = input("\nReactome steps도 등록하시겠습니까? (y/n): ").strip().lower()
    reactome_steps = []

    if add_reactome == "y":
        reactome_steps = _manual_reactome_steps()

    reactome_keywords_text = input(
        "Reactome keyword들을 쉼표로 입력, 기본 target만 사용: "
    ).strip()

    reactome_keywords = [target_text]
    if reactome_keywords_text:
        reactome_keywords.extend(
            [x.strip() for x in reactome_keywords_text.split(",") if x.strip()]
        )

    config = {
        "aliases": aliases,
        "display_name": display_name,
    }

    if kegg_config:
        config["kegg"] = kegg_config

    config["reactome"] = {
        "steps": reactome_steps,
        "keywords": reactome_keywords,
        "representative_paths": [],
    }

    print("\n생성될 registry key:", target_key)
    print("KEGG 등록:", bool(kegg_config))
    print("Reactome step 수:", len(reactome_steps))

    save_answer = input("이 설정을 custom registry에 저장하시겠습니까? (y/n): ").strip().lower()

    if save_answer != "y":
        print("저장하지 않았습니다.")
        return None

    save_custom_target(target_key, config)

    print(f"저장 완료: workflows/target_registry_custom.json")
    print(f"다음 실행부터 target '{target_text}'은/는 자동 workflow로 인식됩니다.")

    return target_key
