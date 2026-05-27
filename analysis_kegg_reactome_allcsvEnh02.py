# -*- coding: utf-8 -*-
"""
analysis_kegg_reactome.py

KEGG/Reactome HIF-1 workflow 결과 ZIP 또는 결과 폴더를 읽어서
다음 6개 관점의 분석 보고서를 생성한다.

1. 전체 규모에서 유추할 수 있는 것
2. KEGG 결과에서 유추할 수 있는 것
3. Reactome 결과에서 유추할 수 있는 것
4. Reactome path 결과에서 유추할 수 있는 것
5. KEGG vs Reactome에서 직관적으로 비교되는 차이
6. HIF-1 생물학 관점에서 바로 유추 가능한 구조

사용 예:
  python analysis_kegg_reactome.py ^
    --kegg-zip results/kegg_hif1_hsa04066_20260525_180659.zip ^
    --reactome-zip results/reactome_hif1_combined_20260525_180701.zip ^
    --summary-zip results/hif1_combined_summary_20260525_181007.zip ^
    --output-dir results/analysis_hif1

또는 이미 압축을 푼 폴더라면:
  python analysis_kegg_reactome.py ^
    --kegg-dir results/kegg_hif1_hsa04066_20260525_180659 ^
    --reactome-dir results/reactome_hif1_combined_20260525_180701 ^
    --summary-dir results/hif1_combined_summary_20260525_181007 ^
    --output-dir results/analysis_hif1
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd


# ---------------------------------------------------------------------
# File loading helpers
# ---------------------------------------------------------------------

def unzip_if_needed(zip_path: str | None, work_dir: Path) -> Path | None:
    if not zip_path:
        return None

    zip_path = Path(zip_path)
    if not zip_path.exists():
        raise FileNotFoundError(f"ZIP not found: {zip_path}")

    extract_dir = work_dir / zip_path.stem
    extract_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)

    # ZIP 내부가 root_folder/file 구조이면 root_folder를 반환
    children = [p for p in extract_dir.iterdir()]
    dirs = [p for p in children if p.is_dir()]
    files = [p for p in children if p.is_file()]

    if len(dirs) == 1 and not files:
        return dirs[0]

    return extract_dir


def resolve_input_dir(
    dir_path: str | None,
    zip_path: str | None,
    work_dir: Path,
) -> Path | None:
    if dir_path:
        p = Path(dir_path)
        if not p.exists():
            raise FileNotFoundError(f"Directory not found: {p}")
        return p

    return unzip_if_needed(zip_path, work_dir)


def read_json(path: str | Path) -> dict[str, Any] | list[Any]:
    """
    JSON 파일을 읽어 dict 또는 list로 반환한다.

    반드시 존재해야 하는 workflow_summary.json, combined_workflow_summary.json
    등을 읽을 때 사용한다.
    """
    json_path = Path(path)

    if not json_path.exists():
        raise FileNotFoundError(f"JSON file not found: {json_path}")

    if not json_path.is_file():
        raise FileNotFoundError(f"Not a file: {json_path}")

    try:
        with json_path.open("r", encoding="utf-8") as f:
            return json.load(f)

    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON file: {json_path}") from exc


def read_json_if_exists(
    path: str | Path,
    default: Any | None = None,
) -> Any:
    """
    JSON 파일이 있으면 읽고, 없으면 default를 반환한다.
    선택적 summary 파일을 읽을 때 사용한다.
    """
    json_path = Path(path)

    if default is None:
        default = {}

    if not json_path.exists() or not json_path.is_file():
        return default

    try:
        return read_json(json_path)
    except Exception as exc:
        print(f"[warning] JSON read failed: {json_path} ({exc})")
        return default


def read_csv_optional(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()

    try:
        return pd.read_csv(path, encoding="utf-8-sig")
    except UnicodeDecodeError:
        return pd.read_csv(path)


def read_json_optional(path: Path) -> dict[str, Any]:
    data = read_json_if_exists(path, default={})
    return data if isinstance(data, dict) else {}


def find_latest_result_dir(
    results_root: Path,
    pattern: str,
) -> Path | None:
    """
    results/ 아래에서 pattern에 맞는 가장 최근 수정 디렉토리를 찾는다.
    """
    if not results_root.exists():
        return None

    candidates = [path for path in results_root.glob(pattern) if path.is_dir()]

    if not candidates:
        return None

    return max(candidates, key=lambda p: p.stat().st_mtime)


def resolve_or_find_input_dir(
    dir_path: str | None,
    zip_path: str | None,
    work_dir: Path,
    *,
    results_root: Path,
    auto_pattern: str,
    label: str,
) -> Path | None:
    """
    1. dir_path 또는 zip_path가 있으면 기존 resolve_input_dir() 사용
    2. 둘 다 없으면 results_root에서 auto_pattern으로 최신 디렉토리 검색
    """
    resolved = resolve_input_dir(
        dir_path=dir_path,
        zip_path=zip_path,
        work_dir=work_dir,
    )

    if resolved is not None:
        return resolved

    latest = find_latest_result_dir(results_root, auto_pattern)

    if latest is not None:
        print(f"[auto] {label} result directory found: {latest}")
        return latest

    print(f"[warning] {label} result directory not found. pattern={auto_pattern}")
    return None


def read_all_csv_files(input_dir: Path) -> dict[str, pd.DataFrame]:
    """
    input_dir 내부의 모든 csv 파일을 읽어서 {파일명: DataFrame} 형태로 반환한다.
    """
    csv_tables: dict[str, pd.DataFrame] = {}

    if input_dir is None or not input_dir.exists():
        return csv_tables

    for csv_path in sorted(input_dir.glob("*.csv")):
        try:
            try:
                csv_tables[csv_path.name] = pd.read_csv(csv_path, encoding="utf-8-sig")
            except UnicodeDecodeError:
                csv_tables[csv_path.name] = pd.read_csv(csv_path)
        except Exception as exc:
            print(f"[warning] CSV read failed: {csv_path} ({exc})")

    return csv_tables


def _find_first_matching_table(
    csv_tables: dict[str, pd.DataFrame],
    patterns: list[str],
) -> pd.DataFrame:
    """
    파일명에 patterns가 모두 포함되는 첫 번째 CSV table을 반환한다.
    """
    for filename, df in csv_tables.items():
        lower = filename.lower()
        if all(pattern.lower() in lower for pattern in patterns):
            return df

    return pd.DataFrame()


def value_counts_dict(df: pd.DataFrame, column: str) -> dict[str, int]:
    if df.empty or column not in df.columns:
        return {}

    return {
        str(k): int(v)
        for k, v in df[column].fillna("").replace("", "EMPTY").value_counts().items()
    }


def contains_text(series: pd.Series, keyword: str) -> pd.Series:
    return series.fillna("").astype(str).str.contains(keyword, case=False, regex=False)


def df_to_records(df: pd.DataFrame, max_rows: int = 20) -> list[dict[str, Any]]:
    if df.empty:
        return []

    return df.head(max_rows).fillna("").to_dict(orient="records")


def unique_path_count(df: pd.DataFrame) -> int:
    """path_index 기준 실제 path 개수를 계산한다."""
    if df.empty or "path_index" not in df.columns:
        return 0
    return int(df["path_index"].nunique())


def path_row_count(df: pd.DataFrame) -> int:
    """path CSV의 step row 개수를 반환한다."""
    return int(len(df))


def nonempty_count(df: pd.DataFrame, column: str) -> int:
    if df.empty or column not in df.columns:
        return 0
    return int(df[column].fillna("").astype(str).str.strip().ne("").sum())


def table_filter_contains(df: pd.DataFrame, keyword: str, columns: list[str] | None = None) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    if columns is None:
        columns = list(df.columns)
    columns = [c for c in columns if c in df.columns]
    if not columns:
        return pd.DataFrame()
    mask = pd.Series(False, index=df.index)
    for col in columns:
        mask = mask | contains_text(df[col], keyword)
    return df[mask]


def build_quality_flags(
    kegg_analysis: dict[str, Any],
    reactome_path_analysis: dict[str, Any],
) -> dict[str, Any]:
    """보고서에서 주의해야 할 사항을 자동 flag로 정리한다."""
    return {
        "kegg_group_expansion_detected": bool(kegg_analysis.get("group_count", 0) > 0),
        "kegg_via_group_edges_present": bool(kegg_analysis.get("via_group_edge_count", 0) > 0),
        "kegg_directed_hif1a_vegfa_exists": bool(kegg_analysis.get("directed_HIF1A_to_VEGFA_path_count", 0) > 0),
        "kegg_directed_hif1a_vegfa_has_via_group": bool(kegg_analysis.get("directed_HIF1A_to_VEGFA_via_group", False)),
        "kegg_compound_name_enrichment_complete": bool(
            kegg_analysis.get("compound_count", 0) > 0
            and kegg_analysis.get("compound_name_missing_count", 0) == 0
        ),
        "kegg_compound_name_enrichment_partial_or_missing": bool(
            kegg_analysis.get("compound_count", 0) > 0
            and kegg_analysis.get("compound_name_missing_count", 0) > 0
        ),
        "kegg_undirected_path_too_many": bool(kegg_analysis.get("undirected_HIF1A_to_VEGFA_path_count", 0) > 100),
        "reactome_directed_path_absent": bool(reactome_path_analysis.get("directed_path_count", 0) == 0),
        "reactome_undirected_path_too_many": bool(reactome_path_analysis.get("undirected_path_count", 0) > 100),
    }


def build_result_file_guide(loaded_csv_files: dict[str, list[str]]) -> dict[str, list[str]]:
    """읽어온 CSV 파일 목록을 기반으로 어떤 파일을 보면 좋은지 안내한다."""
    kegg_files = loaded_csv_files.get("kegg", [])
    reactome_files = loaded_csv_files.get("reactome", [])

    kegg_guide: list[str] = []
    reactome_guide: list[str] = []

    def add_if_present(files: list[str], guide: list[str], filename: str, desc: str) -> None:
        if filename in files:
            guide.append(f"{filename}: {desc}")

    add_if_present(kegg_files, kegg_guide, "kegg_hif1_keyword_hits.csv", "HIF1A, VEGFA, MTOR, ARNT, VHL 등 주요 KEGG node 확인")
    add_if_present(kegg_files, kegg_guide, "activation_edges.csv", "activation 관계 확인")
    add_if_present(kegg_files, kegg_guide, "inhibition_edges.csv", "inhibition 관계 확인")
    add_if_present(kegg_files, kegg_guide, "expression_edges.csv", "expression 및 via_group expression 관계 확인")
    add_if_present(kegg_files, kegg_guide, "kegg_compounds.csv", "KEGG compound ID/name 보강 상태 확인")
    add_if_present(kegg_files, kegg_guide, "paths_directed_HIF1A_to_VEGFA.csv", "HIF1A→VEGFA directed path 및 via_group 여부 확인")
    add_if_present(kegg_files, kegg_guide, "paths_undirected_HIF1A_to_VEGFA.csv", "HIF1A–VEGFA 전체 연결성 확인, 과도한 simple path 해석 주의")

    add_if_present(reactome_files, reactome_guide, "reactome_hif1_keyword_hits.csv", "HIF1A, VEGFA, ARNT, O2, 2OG, CO2, VHL, PHD 확인")
    add_if_present(reactome_files, reactome_guide, "reaction_io.csv", "reaction input/output 해석의 핵심")
    add_if_present(reactome_files, reactome_guide, "event_order_links.csv", "event 순서, precedingEvent 확인")
    add_if_present(reactome_files, reactome_guide, "pathway_event_links.csv", "Pathway가 포함하는 Event 확인")
    add_if_present(reactome_files, reactome_guide, "chemical_nodes.csv", "HIF oxygen-sensing 관련 chemical node 확인")
    add_if_present(reactome_files, reactome_guide, "complex_nodes.csv", "HIF, VHL complex 등 complex node 확인")
    add_if_present(reactome_files, reactome_guide, "entity_pathway_links.csv", "HIF1A가 등장하는 pathway context 확인")

    return {"kegg": kegg_guide, "reactome": reactome_guide}


# ---------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------

def load_kegg_results(kegg_dir: Path) -> dict[str, Any]:
    """
    KEGG 결과 디렉토리 안의 모든 CSV를 먼저 읽고,
    분석에 필요한 대표 table을 파일명 pattern으로 연결한다.
    """
    csv_tables = read_all_csv_files(kegg_dir)

    return {
        "dir": str(kegg_dir),
        "csv_tables": csv_tables,
        "nodes": csv_tables.get("nodes.csv", pd.DataFrame()),
        "edges": csv_tables.get("edges.csv", pd.DataFrame()),
        "keyword_hits": _find_first_matching_table(
            csv_tables,
            patterns=["keyword_hits"],
        ),
        "compounds": _find_first_matching_table(
            csv_tables,
            patterns=["compound"],
        ),
        "directed_paths": _find_first_matching_table(
            csv_tables,
            patterns=["paths_directed"],
        ),
        "undirected_paths": _find_first_matching_table(
            csv_tables,
            patterns=["paths_undirected"],
        ),
        "activation_edges": _find_first_matching_table(
            csv_tables,
            patterns=["activation"],
        ),
        "inhibition_edges": _find_first_matching_table(
            csv_tables,
            patterns=["inhibition"],
        ),
        "expression_edges": _find_first_matching_table(
            csv_tables,
            patterns=["expression"],
        ),
        "indirect_edges": _find_first_matching_table(
            csv_tables,
            patterns=["indirect"],
        ),
        "summary": read_json_if_exists(kegg_dir / "workflow_summary.json", default={}),
    }


def load_reactome_results(reactome_dir: Path) -> dict[str, Any]:
    """
    Reactome 결과 디렉토리 안의 모든 CSV를 먼저 읽고,
    분석에 필요한 대표 table을 파일명 pattern으로 연결한다.
    """
    csv_tables = read_all_csv_files(reactome_dir)

    return {
        "dir": str(reactome_dir),
        "csv_tables": csv_tables,
        "nodes": csv_tables.get("nodes.csv", pd.DataFrame()),
        "edges": csv_tables.get("edges.csv", pd.DataFrame()),
        "keyword_hits": _find_first_matching_table(
            csv_tables,
            patterns=["keyword_hits"],
        ),
        "chemical_nodes": csv_tables.get("chemical_nodes.csv", pd.DataFrame()),
        "complex_nodes": csv_tables.get("complex_nodes.csv", pd.DataFrame()),
        "reaction_event_nodes": csv_tables.get("reaction_event_nodes.csv", pd.DataFrame()),
        "pathway_nodes": csv_tables.get("pathway_nodes.csv", pd.DataFrame()),
        "pathway_event_links": csv_tables.get("pathway_event_links.csv", pd.DataFrame()),
        "reaction_io": csv_tables.get("reaction_io.csv", pd.DataFrame()),
        "catalyst_regulator_links": csv_tables.get("catalyst_regulator_links.csv", pd.DataFrame()),
        "entity_complex_links": csv_tables.get("entity_complex_links.csv", pd.DataFrame()),
        "entity_pathway_links": csv_tables.get("entity_pathway_links.csv", pd.DataFrame()),
        "entity_event_links": csv_tables.get("entity_event_links.csv", pd.DataFrame()),
        "event_order_links": csv_tables.get("event_order_links.csv", pd.DataFrame()),
        "directed_paths": _find_first_matching_table(
            csv_tables,
            patterns=["paths_directed"],
        ),
        "undirected_paths": _find_first_matching_table(
            csv_tables,
            patterns=["paths_undirected"],
        ),
        "summary": read_json_if_exists(reactome_dir / "workflow_summary.json", default={}),
    }

def analyze_overall(kegg: dict[str, Any], reactome: dict[str, Any]) -> dict[str, Any]:
    k_nodes = kegg["nodes"]
    k_edges = kegg["edges"]
    r_nodes = reactome["nodes"]
    r_edges = reactome["edges"]

    return {
        "kegg_node_count": int(len(k_nodes)),
        "kegg_edge_count": int(len(k_edges)),
        "reactome_node_count": int(len(r_nodes)),
        "reactome_edge_count": int(len(r_edges)),
        "kegg_node_types": value_counts_dict(k_nodes, "type"),
        "kegg_relation_types": value_counts_dict(k_edges, "relation_type"),
        "reactome_node_types": value_counts_dict(r_nodes, "type"),
        "reactome_relation_types": value_counts_dict(r_edges, "relation_type"),
        "inference": [
            "KEGG는 gene/compound/map 중심의 signaling overview에 적합하다.",
            "Reactome은 Reaction/Event/Entity/Complex/Chemical 상태 그래프에 적합하다.",
            "Reactome은 node 수가 더 많아 HIF1A/VEGFA/chemical/complex의 상태와 event를 더 세분화해 표현한다. 다만 KEGG edge 수는 group expansion과 via_group edge 추가로 증가할 수 있으므로, edge 수만으로 생물학적 복잡도를 비교하면 안 된다.",
        ],
    }


def analyze_kegg(kegg: dict[str, Any]) -> dict[str, Any]:
    nodes = kegg["nodes"]
    edges = kegg["edges"]
    keyword_hits = kegg["keyword_hits"]
    compounds = kegg["compounds"]
    directed_paths = kegg["directed_paths"]
    undirected_paths = kegg["undirected_paths"]

    def hit_count(keyword: str) -> int:
        if keyword_hits.empty or "keyword" not in keyword_hits.columns:
            return 0
        return int((keyword_hits["keyword"].fillna("").astype(str).str.upper() == keyword.upper()).sum())

    core_keywords = [
        "HIF1A", "VEGF", "VEGFA", "MTOR", "ARNT", "VHL", "EGLN", "PHD", "PDK1", "SLC2A1"
    ]

    hif1a_exists = False
    vegfa_exists = False

    if "node_id" in nodes.columns:
        hif1a_exists = bool((nodes["node_id"].astype(str) == "KEGG:3091").any())
        vegfa_exists = bool((nodes["node_id"].astype(str) == "KEGG:7422").any())

    core_hit_counts = {kw: hit_count(kw) for kw in core_keywords}

    hif1a_edges = pd.DataFrame()
    if not edges.empty and {"source_node", "target_node"}.issubset(edges.columns):
        hif1a_edges = edges[
            (edges["source_node"].astype(str) == "KEGG:3091")
            | (edges["target_node"].astype(str) == "KEGG:3091")
        ]

    group_nodes = pd.DataFrame()
    if not nodes.empty and "type" in nodes.columns:
        group_nodes = nodes[nodes["type"].fillna("").astype(str).str.lower() == "group"]

    via_group_edges = pd.DataFrame()
    if not edges.empty and "relation_type" in edges.columns:
        via_group_edges = edges[contains_text(edges["relation_type"], "via_group")]

    group_component_edges = pd.DataFrame()
    if not edges.empty and "relation_type" in edges.columns:
        group_component_edges = edges[edges["relation_type"].fillna("").astype(str) == "group_component"]

    # group_components/group_component_count 컬럼이 nodes.csv에 없을 경우,
    # group_component edge를 이용해 보고서용 컬럼을 보강한다.
    if not group_nodes.empty:
        group_nodes = group_nodes.copy()
        if "group_components" not in group_nodes.columns:
            group_nodes["group_components"] = ""
        if "group_component_count" not in group_nodes.columns:
            group_nodes["group_component_count"] = 0

        if not group_component_edges.empty and {"source_node", "target_node"}.issubset(group_component_edges.columns):
            for idx, row in group_nodes.iterrows():
                group_id = str(row.get("node_id", ""))
                component_ids = group_component_edges.loc[
                    group_component_edges["source_node"].astype(str) == group_id,
                    "target_node",
                ].fillna("").astype(str).tolist()
                if component_ids:
                    group_nodes.at[idx, "group_components"] = ";".join(component_ids)
                    group_nodes.at[idx, "group_component_count"] = len(component_ids)

    # compound_id/compound_name 컬럼이 kegg_compounds.csv에 없을 경우,
    # raw_id와 label을 이용해 보고서용 컬럼을 보강한다.
    if not compounds.empty:
        compounds = compounds.copy()
        if "compound_id" not in compounds.columns:
            if "raw_id" in compounds.columns:
                compounds["compound_id"] = compounds["raw_id"].fillna("").astype(str)
            elif "node_id" in compounds.columns:
                compounds["compound_id"] = compounds["node_id"].fillna("").astype(str).str.replace("KEGG:", "", regex=False)
            else:
                compounds["compound_id"] = ""
        if "compound_name" not in compounds.columns:
            if "label" in compounds.columns:
                labels = compounds["label"].fillna("").astype(str).str.strip()
                compounds["compound_name"] = labels.where(~labels.str.match(r"^C\d{5}$"), "")
            else:
                compounds["compound_name"] = ""

    directed_via_group = False
    if not directed_paths.empty:
        if "relation_to_next" in directed_paths.columns:
            directed_via_group = bool(contains_text(directed_paths["relation_to_next"], "via_group").any())
        elif "relation_type" in directed_paths.columns:
            directed_via_group = bool(contains_text(directed_paths["relation_type"], "via_group").any())
        elif not edges.empty and "relation_type" in edges.columns and {"source_node", "target_node"}.issubset(edges.columns):
            direct_edge = edges[
                (edges["source_node"].astype(str) == "KEGG:3091")
                & (edges["target_node"].astype(str) == "KEGG:7422")
            ]
            directed_via_group = bool((not direct_edge.empty) and contains_text(direct_edge["relation_type"], "via_group").any())

    compound_name_filled_count = 0
    compound_name_missing_count = int(len(compounds))
    if not compounds.empty:
        if "compound_name" in compounds.columns:
            compound_name_filled_count = nonempty_count(compounds, "compound_name")
            compound_name_missing_count = int(len(compounds) - compound_name_filled_count)
        elif "label" in compounds.columns:
            # label이 Cxxxxx 그대로가 아닌 경우를 보강 성공으로 간주하는 fallback
            labels = compounds["label"].fillna("").astype(str).str.strip()
            compound_name_filled_count = int((labels != "").sum() - labels.str.match(r"^C\d{5}$").sum())
            compound_name_missing_count = int(len(compounds) - max(compound_name_filled_count, 0))

    if path_row_count(directed_paths) > 0:
        if directed_via_group:
            directed_inference = (
                "현재 KEGG 결과에서는 HIF1A→VEGFA directed path가 존재하지만, "
                "이는 group relation을 component node로 확장한 via_group 기반 inferred path로 해석해야 한다."
            )
        else:
            directed_inference = "현재 KEGG 결과에서는 HIF1A→VEGFA directed path가 존재한다."
    else:
        directed_inference = (
            "현재 KEGG directed graph에서는 HIF1A→VEGFA path가 비어 있으며, "
            "이는 KGML relation 방향상 해당 chain이 명시되지 않았다는 뜻이다."
        )

    if compound_name_missing_count == 0 and len(compounds) > 0:
        compound_inference = "KEGG compound ID가 compound_name으로 보강되어 생화학적 해석이 개선되었다."
    else:
        compound_inference = "KEGG compound node 중 compound_name이 비어 있는 항목이 있어 추가 보강 여부를 확인해야 한다."

    return {
        "node_count": int(len(nodes)),
        "edge_count": int(len(edges)),
        "node_types": value_counts_dict(nodes, "type"),
        "relation_types": value_counts_dict(edges, "relation_type"),
        "core_hit_counts": core_hit_counts,
        "HIF1A_node_exists": hif1a_exists,
        "VEGFA_node_exists": vegfa_exists,
        "compound_count": int(len(compounds)),
        "compound_name_filled_count": int(compound_name_filled_count),
        "compound_name_missing_count": int(compound_name_missing_count),
        "compound_samples": df_to_records(compounds, 20),
        "group_count": int(len(group_nodes)),
        "group_samples": df_to_records(group_nodes, 20),
        "group_component_edge_count": int(len(group_component_edges)),
        "via_group_edge_count": int(len(via_group_edges)),
        "via_group_edge_samples": df_to_records(via_group_edges, 20),
        "directed_HIF1A_to_VEGFA_row_count": path_row_count(directed_paths),
        "directed_HIF1A_to_VEGFA_path_count": unique_path_count(directed_paths),
        "directed_HIF1A_to_VEGFA_via_group": directed_via_group,
        "directed_HIF1A_to_VEGFA_path_sample": df_to_records(directed_paths, 20),
        "undirected_HIF1A_to_VEGFA_row_count": path_row_count(undirected_paths),
        "undirected_HIF1A_to_VEGFA_path_count": unique_path_count(undirected_paths),
        "HIF1A_neighbor_edges_sample": df_to_records(hif1a_edges, 30),
        "inference": [
            directed_inference,
            "via_group edge는 원본 KGML direct relation이 아니라 group component 확장으로 생성된 inferred edge이므로 직접 조절 관계와 구분해야 한다.",
            "group_component edge와 group node 개수는 KGML group entry가 실제 component node로 풀렸는지 확인하는 지표이다.",
            "undirected path row 수는 step 단위 row 수이며, 실제 path 수는 unique path_index 기준으로 해석해야 한다.",
            compound_inference,
        ],
    }

def analyze_reactome(reactome: dict[str, Any]) -> dict[str, Any]:
    nodes = reactome["nodes"]
    edges = reactome["edges"]
    keyword_hits = reactome["keyword_hits"]
    chemicals = reactome["chemical_nodes"]
    complexes = reactome["complex_nodes"]
    reaction_io = reactome["reaction_io"]
    entity_pathway_links = reactome["entity_pathway_links"]
    pathway_event_links = reactome["pathway_event_links"]
    event_order_links = reactome["event_order_links"]

    def node_hits(keyword: str) -> list[dict[str, Any]]:
        if nodes.empty:
            return []
        search_cols = [c for c in ["node_id", "label", "raw_id", "reference_identifier", "reference_display"] if c in nodes.columns]
        if not search_cols:
            return []
        mask = pd.Series(False, index=nodes.index)
        for c in search_cols:
            mask = mask | contains_text(nodes[c], keyword)
        return df_to_records(nodes[mask], 20)

    hif1a_hits = node_hits("HIF1A")
    vegfa_hits = node_hits("VEGFA")
    arnt_hits = node_hits("ARNT")

    vegfa_related = pd.DataFrame()
    if not reaction_io.empty:
        cols = [c for c in reaction_io.columns if "label" in c or "node" in c]
        mask = pd.Series(False, index=reaction_io.index)
        for c in cols:
            mask = mask | contains_text(reaction_io[c], "VEGFA")
        vegfa_related = reaction_io[mask]

    hif1a_pathway_context = pd.DataFrame()
    if not entity_pathway_links.empty:
        cols = [c for c in entity_pathway_links.columns if "label" in c or "node" in c]
        mask = pd.Series(False, index=entity_pathway_links.index)
        for c in cols:
            mask = mask | contains_text(entity_pathway_links[c], "HIF1A")
        hif1a_pathway_context = entity_pathway_links[mask]

    # R-HSA-1234171: HIF-alpha binds ARNT reaction의 input/output을 명시적으로 추출
    hif_arnt_io = pd.DataFrame()
    if not reaction_io.empty:
        cols = [c for c in reaction_io.columns if "label" in c or "node" in c]
        mask = pd.Series(False, index=reaction_io.index)
        for c in cols:
            mask = mask | contains_text(reaction_io[c], "R-HSA-1234171")
            mask = mask | contains_text(reaction_io[c], "HIF-alpha binds ARNT")
        hif_arnt_io = reaction_io[mask]

    return {
        "node_count": int(len(nodes)),
        "edge_count": int(len(edges)),
        "node_types": value_counts_dict(nodes, "type"),
        "relation_types": value_counts_dict(edges, "relation_type"),
        "chemical_count": int(len(chemicals)),
        "complex_count": int(len(complexes)),
        "reaction_io_count": int(len(reaction_io)),
        "pathway_event_link_count": int(len(pathway_event_links)),
        "event_order_link_count": int(len(event_order_links)),
        "HIF1A_node_samples": hif1a_hits,
        "VEGFA_node_samples": vegfa_hits,
        "ARNT_node_samples": arnt_hits,
        "VEGFA_reaction_io_samples": df_to_records(vegfa_related, 30),
        "HIF_ARNT_reaction_io_samples": df_to_records(hif_arnt_io, 30),
        "HIF1A_pathway_context_samples": df_to_records(hif1a_pathway_context, 30),
        "chemical_samples": df_to_records(chemicals, 30),
        "inference": [
            "Reactome은 HIF1A를 cytosol/nucleoplasm/modified form 등 상태별 entity로 분리한다.",
            "Reactome은 VEGFA gene, Expression of VEGFA event, VEGFA protein product를 분리해 보여줄 수 있다.",
            "HIF-alpha와 ARNT 결합 반응은 reaction input/output으로 해석하는 것이 적절하다.",
            "O2, 2OG, CO2, succinate 같은 chemical node가 SimpleEntity로 나타나 HIF oxygen-sensing mechanism 해석에 유리하다.",
            "entity_pathway_links는 HIF1A가 hypoxia 외의 여러 pathway context에도 등장하는지 확인하는 데 유용하다.",
        ],
    }

def analyze_reactome_paths(reactome: dict[str, Any]) -> dict[str, Any]:
    directed = reactome["directed_paths"]
    undirected = reactome["undirected_paths"]

    def path_count(df: pd.DataFrame) -> int:
        if df.empty or "path_index" not in df.columns:
            return 0
        return int(df["path_index"].nunique())

    relation_counts = {}
    if not undirected.empty and "relation_to_next" in undirected.columns:
        relation_counts = value_counts_dict(undirected, "relation_to_next")

    return {
        "directed_path_row_count": int(len(directed)),
        "directed_path_count": path_count(directed),
        "undirected_path_row_count": int(len(undirected)),
        "undirected_path_count": path_count(undirected),
        "undirected_relation_to_next_counts": relation_counts,
        "undirected_path_sample": df_to_records(undirected, 40),
        "inference": [
            "Reactome directed path가 비어 있으면 엄격한 reaction input/output 방향으로는 직접 흐름이 없다는 뜻이다.",
            "Reactome undirected path가 매우 많으면 연결성은 풍부하지만 생물학적 causality 해석에는 과도할 수 있다.",
            "relation_to_next의 EMPTY는 보통 path의 마지막 node이거나 relation_to_next를 계산할 수 없는 step을 의미하므로 별도 생물학적 relation으로 해석하면 안 된다.",
            "Reactome path 해석은 all simple paths보다 relation_type 필터링, event_order_links, reaction_io 중심으로 보는 것이 좋다.",
        ],
    }


def compare_kegg_reactome(kegg_analysis: dict[str, Any], reactome_analysis: dict[str, Any]) -> dict[str, Any]:
    kegg_chemical = "compound ID/name 보강 상태를 확인해야 함"
    if kegg_analysis.get("compound_count", 0) > 0 and kegg_analysis.get("compound_name_missing_count", 0) == 0:
        kegg_chemical = "compound ID가 name으로 보강됨"

    kegg_path = "directed/undirected signaling relation"
    if kegg_analysis.get("directed_HIF1A_to_VEGFA_via_group"):
        kegg_path = "directed path 존재, 단 via_group inferred edge 해석 주의"

    return {
        "table": [
            {
                "item": "전체 성격",
                "KEGG": "signaling map / gene-compound-relation overview",
                "Reactome": "event-reaction-entity-state graph",
            },
            {
                "item": "주요 node",
                "KEGG": "gene, compound, map, group",
                "Reactome": "Pathway, Reaction, BlackBoxEvent, Protein, Complex, SimpleEntity",
            },
            {
                "item": "HIF1A 표현",
                "KEGG": "대체로 KEGG:3091 하나, group 확장 시 downstream inferred edge 증가",
                "Reactome": "cytosol/nucleoplasm/modified form으로 분리 가능",
            },
            {
                "item": "VEGFA 표현",
                "KEGG": "gene node 및 via_group expression edge로 연결 가능",
                "Reactome": "VEGFA gene → Expression event → VEGFA product",
            },
            {
                "item": "Chemical 해석",
                "KEGG": kegg_chemical,
                "Reactome": "ChEBI/reference/compartment까지 해석 가능",
            },
            {
                "item": "path 해석",
                "KEGG": kegg_path,
                "Reactome": "reaction I/O, precedingEvent, pathway-event 관계 중심",
            },
        ],
        "inference": [
            "KEGG는 빠른 pathway overview에 적합하고 Reactome은 mechanism/state 해석에 적합하다.",
            "KEGG의 via_group edge는 group expansion 기반 inferred relation이므로 Reactome의 reaction I/O와 같은 수준의 직접 기전으로 해석하면 안 된다.",
            "두 결과는 상호 대체가 아니라 상호 보완적으로 보는 것이 좋다.",
        ],
    }

def infer_hif1_biology(kegg_analysis: dict[str, Any], reactome_analysis: dict[str, Any]) -> dict[str, Any]:
    return {
        "candidate_structure": [
            "저산소 또는 산소 반응 context",
            "O2, 2OG 의존적 PHD hydroxylation",
            "HIF1A/HIF2A/HIF3A의 compartment별 변형",
            "VHL/EloB,C/CUL2/RBX1 관련 ubiquitination/degradation",
            "HIF-alpha nuclear translocation",
            "HIF-alpha + ARNT 결합",
            "HIF complex 형성",
            "HIF:CBP:p300 promoter complex",
            "VEGFA, EPO, CA9, HIGD1A expression",
        ],
        "inference": [
            "KEGG 결과는 HIF-1 signaling pathway의 전체 구성요소와 relation overview를 준다.",
            "Reactome 결과는 HIF oxygen-sensing, hydroxylation, complex formation, expression event를 상태 중심으로 구체화한다.",
            "HIF1A→VEGFA 생물학적 관계는 directed simple path 하나로만 판단하면 부족하며, expression event와 pathway context를 함께 봐야 한다.",
        ],
    }


# ---------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------

def md_dict_counts(title: str, counts: dict[str, int]) -> str:
    lines = [f"### {title}", ""]
    if not counts:
        lines.append("- 없음")
        return "\n".join(lines)

    lines.append("| 항목 | 개수 |")
    lines.append("|---|---:|")
    for k, v in counts.items():
        lines.append(f"| {k} | {v} |")
    return "\n".join(lines)


def render_records(records: list[dict[str, Any]], max_rows: int = 10) -> str:
    if not records:
        return "- 없음"

    records = records[:max_rows]
    keys = list(records[0].keys())
    lines = []
    lines.append("| " + " | ".join(keys) + " |")
    lines.append("|" + "|".join(["---"] * len(keys)) + "|")

    for rec in records:
        vals = [str(rec.get(k, "")).replace("\n", " ")[:120] for k in keys]
        lines.append("| " + " | ".join(vals) + " |")

    return "\n".join(lines)


def build_markdown_report(
    overall: dict[str, Any],
    kegg_analysis: dict[str, Any],
    reactome_analysis: dict[str, Any],
    reactome_path_analysis: dict[str, Any],
    comparison: dict[str, Any],
    biology: dict[str, Any],
    quality_flags: dict[str, Any] | None = None,
    result_file_guide: dict[str, list[str]] | None = None,
) -> str:
    lines: list[str] = []
    quality_flags = quality_flags or {}
    result_file_guide = result_file_guide or {"kegg": [], "reactome": []}

    lines.append("# HIF-1 KEGG / Reactome 결과 분석 보고서")
    lines.append("")

    # 1
    lines.append("## 1. 전체 규모에서 유추할 수 있는 것")
    lines.append("")
    lines.append(f"- KEGG: nodes={overall['kegg_node_count']}, edges={overall['kegg_edge_count']}")
    lines.append(f"- Reactome: nodes={overall['reactome_node_count']}, edges={overall['reactome_edge_count']}")
    lines.append("")
    lines.append(md_dict_counts("KEGG node type", overall["kegg_node_types"]))
    lines.append("")
    lines.append(md_dict_counts("KEGG relation type", overall["kegg_relation_types"]))
    lines.append("")
    lines.append(md_dict_counts("Reactome node type", overall["reactome_node_types"]))
    lines.append("")
    lines.append(md_dict_counts("Reactome relation type", overall["reactome_relation_types"]))
    lines.append("")
    for x in overall["inference"]:
        lines.append(f"- {x}")
    if quality_flags.get("kegg_group_expansion_detected"):
        lines.append("- KEGG node/edge 수 증가는 group_component 및 via_group edge 확장 결과일 수 있으므로, 원본 KGML relation 수와 구분해야 한다.")
    lines.append("")

    # 2
    lines.append("## 2. KEGG 결과에서 유추할 수 있는 것")
    lines.append("")
    lines.append(f"- HIF1A node exists: {kegg_analysis['HIF1A_node_exists']}")
    lines.append(f"- VEGFA node exists: {kegg_analysis['VEGFA_node_exists']}")
    lines.append(f"- compound count: {kegg_analysis['compound_count']}")
    lines.append(f"- compound name filled/missing: {kegg_analysis.get('compound_name_filled_count', 0)} / {kegg_analysis.get('compound_name_missing_count', 0)}")
    lines.append(f"- group count: {kegg_analysis.get('group_count', 0)}")
    lines.append(f"- group_component edges: {kegg_analysis.get('group_component_edge_count', 0)}")
    lines.append(f"- via_group edges: {kegg_analysis.get('via_group_edge_count', 0)}")
    lines.append(f"- directed HIF1A→VEGFA path count / rows: {kegg_analysis.get('directed_HIF1A_to_VEGFA_path_count', 0)} / {kegg_analysis['directed_HIF1A_to_VEGFA_row_count']}")
    lines.append(f"- directed HIF1A→VEGFA uses via_group: {kegg_analysis.get('directed_HIF1A_to_VEGFA_via_group', False)}")
    lines.append(f"- undirected HIF1A—VEGFA path count / rows: {kegg_analysis.get('undirected_HIF1A_to_VEGFA_path_count', 0)} / {kegg_analysis['undirected_HIF1A_to_VEGFA_row_count']}")
    lines.append("")
    lines.append("### KEGG 핵심 keyword hit count")
    lines.append("")
    lines.append("| keyword | hit count |")
    lines.append("|---|---:|")
    for k, v in kegg_analysis["core_hit_counts"].items():
        lines.append(f"| {k} | {v} |")
    lines.append("")
    lines.append("### KEGG directed HIF1A→VEGFA path sample")
    lines.append("")
    lines.append(render_records(kegg_analysis.get("directed_HIF1A_to_VEGFA_path_sample", []), 20))
    lines.append("")
    lines.append("### KEGG compound samples")
    lines.append("")
    lines.append(render_records(kegg_analysis.get("compound_samples", []), 10))
    lines.append("")
    lines.append("### KEGG group samples")
    lines.append("")
    lines.append(render_records(kegg_analysis.get("group_samples", []), 10))
    lines.append("")
    lines.append("### KEGG via_group edge samples")
    lines.append("")
    lines.append(render_records(kegg_analysis.get("via_group_edge_samples", []), 10))
    lines.append("")
    lines.append("### HIF1A 주변 edge sample")
    lines.append("")
    lines.append(render_records(kegg_analysis["HIF1A_neighbor_edges_sample"], 20))
    lines.append("")
    for x in kegg_analysis["inference"]:
        lines.append(f"- {x}")
    lines.append("")

    # 3
    lines.append("## 3. Reactome 결과에서 유추할 수 있는 것")
    lines.append("")
    lines.append(f"- chemical count: {reactome_analysis['chemical_count']}")
    lines.append(f"- complex count: {reactome_analysis['complex_count']}")
    lines.append(f"- reaction_io rows: {reactome_analysis['reaction_io_count']}")
    lines.append(f"- pathway_event links: {reactome_analysis['pathway_event_link_count']}")
    lines.append(f"- event_order links: {reactome_analysis['event_order_link_count']}")
    lines.append("")
    lines.append("### HIF1A node samples")
    lines.append("")
    lines.append(render_records(reactome_analysis["HIF1A_node_samples"], 10))
    lines.append("")
    lines.append("### VEGFA node samples")
    lines.append("")
    lines.append(render_records(reactome_analysis["VEGFA_node_samples"], 10))
    lines.append("")
    lines.append("### VEGFA reaction I/O samples")
    lines.append("")
    lines.append(render_records(reactome_analysis["VEGFA_reaction_io_samples"], 20))
    lines.append("")
    lines.append("### HIF-alpha binds ARNT reaction I/O samples")
    lines.append("")
    lines.append(render_records(reactome_analysis.get("HIF_ARNT_reaction_io_samples", []), 20))
    lines.append("")
    lines.append("### Chemical node samples")
    lines.append("")
    lines.append(render_records(reactome_analysis["chemical_samples"], 20))
    lines.append("")
    for x in reactome_analysis["inference"]:
        lines.append(f"- {x}")
    lines.append("")

    # 4
    lines.append("## 4. Reactome path 결과에서 유추할 수 있는 것")
    lines.append("")
    lines.append(f"- directed path count: {reactome_path_analysis['directed_path_count']}")
    lines.append(f"- directed path rows: {reactome_path_analysis['directed_path_row_count']}")
    lines.append(f"- undirected path count: {reactome_path_analysis['undirected_path_count']}")
    lines.append(f"- undirected path rows: {reactome_path_analysis['undirected_path_row_count']}")
    lines.append("")
    lines.append(md_dict_counts("Undirected path relation_to_next counts", reactome_path_analysis["undirected_relation_to_next_counts"]))
    lines.append("")
    lines.append("### Reactome undirected path sample")
    lines.append("")
    lines.append(render_records(reactome_path_analysis.get("undirected_path_sample", []), 15))
    lines.append("")
    for x in reactome_path_analysis["inference"]:
        lines.append(f"- {x}")
    lines.append("")

    # 5
    lines.append("## 5. KEGG vs Reactome에서 직관적으로 비교되는 차이")
    lines.append("")
    lines.append("| 항목 | KEGG | Reactome |")
    lines.append("|---|---|---|")
    for row in comparison["table"]:
        lines.append(f"| {row['item']} | {row['KEGG']} | {row['Reactome']} |")
    lines.append("")
    for x in comparison["inference"]:
        lines.append(f"- {x}")
    lines.append("")

    # 6
    lines.append("## 6. HIF-1 생물학 관점에서 바로 유추 가능한 구조")
    lines.append("")
    for i, item in enumerate(biology["candidate_structure"], start=1):
        lines.append(f"{i}. {item}")
    lines.append("")
    for x in biology["inference"]:
        lines.append(f"- {x}")
    lines.append("")

    # 7
    lines.append("## 7. 현재 결과에서 주의해야 할 점")
    lines.append("")
    if quality_flags.get("kegg_via_group_edges_present"):
        lines.append("- KEGG의 `via_group` edge는 parser가 KGML group relation을 component node 수준으로 확장한 inferred edge입니다. 원본 KGML의 직접 relation과 구분해야 합니다.")
    if quality_flags.get("kegg_directed_hif1a_vegfa_has_via_group"):
        lines.append("- 현재 HIF1A→VEGFA directed path는 존재하지만, `expression, via_group` 기반이므로 직접 전사 조절 관계로 단정하지 말고 group 확장 결과로 해석해야 합니다.")
    if quality_flags.get("kegg_compound_name_enrichment_partial_or_missing"):
        lines.append("- KEGG compound 중 compound_name이 비어 있는 항목이 있어 KEGG compound enrichment 반영 여부를 확인해야 합니다.")
    if quality_flags.get("kegg_undirected_path_too_many"):
        lines.append("- KEGG undirected path가 매우 많으므로 전체 simple path를 생물학적 flow로 그대로 해석하지 말고 relation_type과 via_group 여부로 필터링해야 합니다.")
    if quality_flags.get("reactome_undirected_path_too_many"):
        lines.append("- Reactome undirected path가 매우 많으므로 연결성 확인 용도로만 사용하고, flow 해석은 reaction_io/event_order/pathway_event 중심으로 해야 합니다.")
    if quality_flags.get("reactome_directed_path_absent"):
        lines.append("- Reactome directed path 부재는 HIF1A와 ARNT 결합 반응이 무관하다는 뜻이 아니라, 엄격한 reaction input/output 방향으로 직접 simple path가 없다는 뜻입니다.")
    lines.append("")

    # 8
    lines.append("## 8. 결과 파일별로 무엇을 보면 좋은가")
    lines.append("")
    lines.append("### KEGG")
    if result_file_guide.get("kegg"):
        for item in result_file_guide["kegg"]:
            lines.append(f"- {item}")
    else:
        lines.append("- 읽어온 KEGG CSV 파일 목록이 없습니다.")
    lines.append("")
    lines.append("### Reactome")
    if result_file_guide.get("reactome"):
        for item in result_file_guide["reactome"]:
            lines.append(f"- {item}")
    else:
        lines.append("- 읽어온 Reactome CSV 파일 목록이 없습니다.")
    lines.append("")

    # 9
    lines.append("## 9. 자동 품질 점검 결과")
    lines.append("")
    if quality_flags:
        flag_descriptions = {
            "kegg_group_expansion_detected": "KEGG group entry가 component node로 확장되었는지 여부",
            "kegg_via_group_edges_present": "via_group inferred edge가 생성되었는지 여부",
            "kegg_directed_hif1a_vegfa_exists": "HIF1A→VEGFA directed simple path 존재 여부",
            "kegg_directed_hif1a_vegfa_has_via_group": "HIF1A→VEGFA path가 via_group 기반인지 여부",
            "kegg_compound_name_enrichment_complete": "KEGG compound name 보강 완료 여부",
            "kegg_compound_name_enrichment_partial_or_missing": "KEGG compound name 일부 누락 여부",
            "kegg_undirected_path_too_many": "KEGG undirected path가 과도하게 많은지 여부",
            "reactome_directed_path_absent": "Reactome directed simple path 부재 여부",
            "reactome_undirected_path_too_many": "Reactome undirected path가 과도하게 많은지 여부",
        }
        lines.append("| flag | value | 해석 |")
        lines.append("|---|---:|---|")
        for flag, value in quality_flags.items():
            lines.append(f"| {flag} | {value} | {flag_descriptions.get(flag, '')} |")
    else:
        lines.append("- quality_flags 정보가 없습니다.")
    lines.append("")

    return "\n".join(lines)

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kegg-zip", default=None)
    parser.add_argument("--reactome-zip", default=None)
    parser.add_argument("--summary-zip", default=None)
    parser.add_argument("--kegg-dir", default=None)
    parser.add_argument("--reactome-dir", default=None)
    parser.add_argument("--summary-dir", default=None)
    parser.add_argument(
        "--results-root",
        default="results",
        help="자동 탐색할 results 루트 디렉토리. 기본값: results",
    )
    parser.add_argument("--output-dir", default="results/analysis_kegg_reactome")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    work_dir = output_dir / "_unzipped"
    work_dir.mkdir(parents=True, exist_ok=True)

    results_root = Path(args.results_root)

    kegg_dir = resolve_or_find_input_dir(
        dir_path=args.kegg_dir,
        zip_path=args.kegg_zip,
        work_dir=work_dir,
        results_root=results_root,
        auto_pattern="kegg_hif1_hsa04066_*",
        label="KEGG HIF-1",
    )

    reactome_dir = resolve_or_find_input_dir(
        dir_path=args.reactome_dir,
        zip_path=args.reactome_zip,
        work_dir=work_dir,
        results_root=results_root,
        auto_pattern="reactome_hif1_combined_*",
        label="Reactome HIF-1",
    )

    summary_dir = resolve_or_find_input_dir(
        dir_path=args.summary_dir,
        zip_path=args.summary_zip,
        work_dir=work_dir,
        results_root=results_root,
        auto_pattern="hif1_combined_summary_*",
        label="Combined summary",
    )

    if kegg_dir is None:
        raise FileNotFoundError(
            "KEGG result directory를 찾지 못했습니다. "
            "--kegg-dir 또는 --kegg-zip을 지정하거나 "
            "results/kegg_hif1_hsa04066_* 폴더를 확인하세요."
        )

    if reactome_dir is None:
        raise FileNotFoundError(
            "Reactome result directory를 찾지 못했습니다. "
            "--reactome-dir 또는 --reactome-zip을 지정하거나 "
            "results/reactome_hif1_combined_* 폴더를 확인하세요."
        )

    print("\n--- Input directories ---")
    print(f"KEGG     : {kegg_dir}")
    print(f"Reactome : {reactome_dir}")
    print(f"Summary  : {summary_dir}")

    kegg = load_kegg_results(kegg_dir)
    reactome = load_reactome_results(reactome_dir)

    combined_summary = {}
    if summary_dir is not None:
        combined_summary = read_json_if_exists(
            summary_dir / "combined_workflow_summary.json",
            default={},
        )

    overall = analyze_overall(kegg, reactome)
    kegg_analysis = analyze_kegg(kegg)
    reactome_analysis = analyze_reactome(reactome)
    reactome_path_analysis = analyze_reactome_paths(reactome)
    comparison = compare_kegg_reactome(kegg_analysis, reactome_analysis)
    biology = infer_hif1_biology(kegg_analysis, reactome_analysis)

    loaded_csv_files = {
        "kegg": sorted(kegg.get("csv_tables", {}).keys()),
        "reactome": sorted(reactome.get("csv_tables", {}).keys()),
    }
    quality_flags = build_quality_flags(kegg_analysis, reactome_path_analysis)
    result_file_guide = build_result_file_guide(loaded_csv_files)

    report = build_markdown_report(
        overall,
        kegg_analysis,
        reactome_analysis,
        reactome_path_analysis,
        comparison,
        biology,
        quality_flags=quality_flags,
        result_file_guide=result_file_guide,
    )

    report_path = output_dir / "analysis_report.md"
    with report_path.open("w", encoding="utf-8") as f:
        f.write(report)

    analysis_json = {
        "kegg_dir": str(kegg_dir),
        "reactome_dir": str(reactome_dir),
        "summary_dir": str(summary_dir) if summary_dir else None,
        "combined_summary": combined_summary,
        "loaded_csv_files": loaded_csv_files,
        "quality_flags": quality_flags,
        "result_file_guide": result_file_guide,
        "overall": overall,
        "kegg_analysis": kegg_analysis,
        "reactome_analysis": reactome_analysis,
        "reactome_path_analysis": reactome_path_analysis,
        "comparison": comparison,
        "biology": biology,
    }

    json_path = output_dir / "analysis_summary.json"
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(analysis_json, f, ensure_ascii=False, indent=2)

    print(f"Analysis report saved: {report_path}")
    print(f"Analysis summary saved: {json_path}")


if __name__ == "__main__":
    main()
