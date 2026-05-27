"""Interactive CLI entry point for S2SignalPathway."""

from __future__ import annotations

import lxml, json, sys, subprocess

from pathlib import Path
from datetime import datetime
import networkx as nx

from SignalPathway.s2pathway import SignalingPathway, PathwaySource
from workflows.target_registry import recognize_target_key, get_target_config


POST_MODELING_CONFIG = {
    "hif1": {
        "target_keywords": ["HIF1A", "VEGFA", "ARNT"],
        "reactome_sbml_ids": [
            "R-HSA-1234174",  # Cellular response to hypoxia
            "R-HSA-1234158",  # Regulation of gene expression by HIF
            "R-HSA-1234171",  # HIF-alpha binds ARNT
        ],
        "sbml_single_check_id": "R-HSA-1234174",
    }
}

def _select_kegg_pathway(pathway: SignalingPathway) -> str:
    mode = input("Choose KEGG mode [1: ID Direct Input, 2: Keyword Search]: ").strip()

    if mode == "2":
        keyword = input("Enter KEGG keyword to search, e.g. MAPK or HIF1: ").strip()
        results = pathway.search_pathway_id(keyword)

        if not results:
            print(f"\n'{keyword}'에 해당하는 KEGG pathway를 찾지 못했습니다.")
            print("  1. 다른 keyword로 다시 검색")
            print("  2. KEGG pathway ID 직접 입력")
            print("  3. 종료")

            retry_choice = input("선택 [1/2/3]: ").strip()

            if retry_choice == "1":
                return _select_kegg_pathway(pathway)

            if retry_choice == "2":
                return input("Enter KEGG pathway ID, e.g. hsa04066 or map04115: ").strip()

            raise ValueError("No KEGG pathways found.")

        print("\n--- KEGG Search Results ---")
        for index, result in enumerate(results):
            print(f"[{index}] {result['id']}: {result['name']}")

        choice_text = input(
            "\nSelect KEGG pathway index or pathway ID, e.g. 0 or map04115 or hsa04115: "
        ).strip()

        if choice_text.isdigit():
            choice = int(choice_text)

            if choice < 0 or choice >= len(results):
                raise IndexError("Invalid selection index.")

            return results[choice]["id"]

        for result in results:
            result_id = str(result["id"]).strip()

            if choice_text.lower() == result_id.lower():
                return result_id

            if choice_text.lower().replace("map", "hsa", 1) == result_id.lower():
                return result_id

            if choice_text.lower().replace("hsa", "map", 1) == result_id.lower():
                return result_id

        if choice_text.lower().startswith(("map", "hsa")):
            return choice_text

        raise ValueError(
            f"Invalid KEGG selection: {choice_text}. "
            "번호 또는 KEGG pathway ID를 입력하세요."
        )

    return input("Enter KEGG pathway ID, e.g. hsa04010: ").strip()

def _select_reactome_pathway() -> str:
    return input("Enter Reactome pathway/reaction ID, e.g. R-HSA-168256 or R-HSA-1234174 or R-HSA-1234158: ").strip()


def _print_node_preview(pathway: SignalingPathway, limit: int = 30) -> None:
    print("\n--- Node Preview ---")

    nodes = list(pathway.graph.nodes(data=True))

    for node_id, attrs in nodes[:limit]:
        print(
            f"[{attrs.get('source')}] {node_id} | "
            f"type={attrs.get('type')} | label={attrs.get('label')}"
        )

    if len(nodes) > limit:
        print(f"... {len(nodes) - limit} more nodes")


def _interactive_label_search(pathway: SignalingPathway) -> None:
    search_confirm = input("\n노드 label로 검색하시겠습니까? (y/n): ").strip().lower()
    if search_confirm != "y":
        return

    while True:
        query = input("검색할 유전자/단백질/반응 이름을 입력하세요. 종료: q > ").strip()

        if query.lower() == "q":
            break

        found_ids = pathway.find_node_by_label(query)

        if not found_ids:
            print(f"'{query}'에 해당하는 노드를 찾을 수 없습니다.")
            continue

        print(f"\n--- Search Results for '{query}' ---")

        for node_id in found_ids:
            attrs = pathway.graph.nodes[node_id]
            print(
                f"{node_id} | "
                f"source={attrs.get('source')} | "
                f"type={attrs.get('type')} | "
                f"label={attrs.get('label')}"
            )


def _interactive_path_search(pathway: SignalingPathway) -> list:
    print("\n--- 경로 탐색 ---")
    print("노드 ID 예:")
    print("  KEGG     : KEGG:3091 또는 3091")
    print("  Reactome : REACTOME:R-HSA-xxxx 또는 R-HSA-xxxx")

    start_node = input("시작 노드 ID 입력: ").strip()
    end_node = input("도착 노드 ID 입력: ").strip()

    if not start_node or not end_node:
        return []

    start_node_norm = pathway.normalize_node_id(start_node)
    end_node_norm = pathway.normalize_node_id(end_node)

    print(f"Normalized start node: {start_node_norm}")
    print(f"Normalized end node  : {end_node_norm}")

    graph = pathway.graph

    if start_node_norm not in graph:
        print(f"시작 노드가 그래프에 없습니다: {start_node_norm}")
        return []

    if end_node_norm not in graph:
        print(f"도착 노드가 그래프에 없습니다: {end_node_norm}")
        return []

    paths = pathway.get_pathways_between_nodes(
        start_node_norm,
        end_node_norm,
        cutoff=10
    )

    if paths:
        print(f"\n방향성 경로가 있습니다. 총 {len(paths)}개의 경로가 발견되었습니다.")

        for index, path in enumerate(paths, start=1):
            path_names = [
                f"{node}({graph.nodes[node].get('label', 'Unknown')})"
                for node in path
            ]
            print(f"[Directed Path {index}] " + " -> ".join(path_names))

        return paths

    print(f"\n방향성 경로는 없습니다: {start_node_norm} -> {end_node_norm}")

    reverse_paths = []
    if nx.has_path(graph, end_node_norm, start_node_norm):
        reverse_paths = list(nx.all_simple_paths(
            graph,
            source=end_node_norm,
            target=start_node_norm,
            cutoff=10
        ))

        print(f"반대 방향 경로는 있습니다: {end_node_norm} -> {start_node_norm}")
        print(f"반대 방향 경로 수: {len(reverse_paths)}")

    undirected_graph = graph.to_undirected()

    if nx.has_path(undirected_graph, start_node_norm, end_node_norm):
        shortest = nx.shortest_path(
            undirected_graph,
            source=start_node_norm,
            target=end_node_norm
        )

        print("\n무방향 연결성은 있습니다.")
        print("무방향 최단 연결:")

        shortest_names = [
            f"{node}({graph.nodes[node].get('label', 'Unknown')})"
            for node in shortest
        ]
        print(" -> ".join(shortest_names))

    else:
        print("무방향으로도 두 노드는 연결되어 있지 않습니다.")

    return []

def _module_name(short_module: str) -> str:
    """
    main.py 위치에서 modeling package를 어떻게 import할지 결정한다.

    우선순위:
      1. modeling.xxx
      2. S2SignalPathway.modeling.xxx
    """
    if Path("modeling").exists():
        return f"modeling.{short_module}"

    if Path("S2SignalPathway").exists():
        return f"S2SignalPathway.modeling.{short_module}"

    return f"modeling.{short_module}"

def _run_command(command: list[str], *, description: str) -> bool:
    """
    외부 python -m 명령을 실행한다.
    실패해도 전체 main workflow를 중단하지 않고 False를 반환한다.
    """
    print(f"\n--- {description} ---")
    print(" ".join(command))

    try:
        subprocess.run(command, check=True)
        return True
    except subprocess.CalledProcessError as exc:
        print(f"[warning] {description} 실패: {exc}")
        return False
    except FileNotFoundError as exc:
        print(f"[warning] 실행 파일을 찾지 못했습니다: {exc}")
        return False
    
def _save_all_results(
    pathway: SignalingPathway,
    selected_id: str,
    result: dict,
    paths: list
) -> None:
    """
    실행 결과를 results/ 하위에 모두 저장한다.
    save_run_results()는 s2pathway.py의 SignalingPathway 클래스에 있어야 한다.
    """
    save_confirm = input("\n실행 결과 전체를 results 디렉토리에 저장하시겠습니까? (y/n): ").strip().lower()

    if save_confirm != "y":
        return

    run_dir = pathway.save_run_results(
        selected_id=selected_id,
        result=result,
        paths=paths,
        output_root="results"
    )

    print(f"Saved all results: {run_dir}")

def normalize_goal_text(goal: str) -> str:
    return goal.strip().lower().replace("_", "-").replace(" ", "")


def recognize_analysis_goal(goal: str) -> str:
    target_key = recognize_target_key(goal)

    if target_key:
        return target_key

    return "manual"
    
def print_goal_workflow_description(goal_key: str) -> None:
    if goal_key == "manual":
        print("\n등록된 자동 workflow가 없습니다.")
        print("기존 수동 모드로 KEGG 또는 Reactome ID를 직접 입력하여 진행합니다.")
        return

    config = get_target_config(goal_key)

    print(f"\n인식된 분석 목표: {config.get('display_name', goal_key)}")

    if "kegg" in config:
        kegg = config["kegg"]
        print("\nKEGG workflow:")
        print(f"  - {kegg.get('pathway_id')} 중심 KGML graph 생성")
        print(f"  - pathway name: {kegg.get('pathway_name', '')}")
        print("  - gene / compound / relation 분석")
        print("  - 대표 path 확인")
        print("  - compound node 추출")

    if "reactome" in config:
        reactome = config["reactome"]
        steps = reactome.get("steps", [])

        print("\nReactome workflow:")

        if steps:
            for step in steps:
                print(f"  - {step['reactome_id']}: {step.get('description', '')}")
            print("  - Pathway/Event/Entity/Reaction/Complex/ChemicalCompound 분리 추출")
        else:
            print("  - Reactome workflow steps가 아직 등록되어 있지 않습니다.")
            
def select_data_source_for_goal() -> str:
    print("\n데이터 소스를 선택하세요.")
    print("  1. KEGG")
    print("  2. Reactome")
    print("  3. Both")
    print("  4. Manual mode")

    choice = input("선택 [1/2/3/4]: ").strip().lower()

    if choice in {"1", "kegg"}:
        return "kegg"

    if choice in {"2", "reactome"}:
        return "reactome"

    if choice in {"3", "both", "all"}:
        return "both"

    return "manual"

def run_goal_workflow(goal_key: str, source_choice: str):
    if goal_key == "manual":
        return None

    from workflows.target_registry import get_target_config
    from workflows.generic_kegg_target_workflow import run_generic_kegg_target_workflow
    from workflows.generic_reactome_target_workflow import run_generic_reactome_target_workflow

    config = get_target_config(goal_key)
    results = {}

    if source_choice == "kegg":
        if "kegg" not in config:
            print(f"{goal_key}에 대한 KEGG workflow 설정이 없습니다.")
            return None

        print(f"\n--- Running KEGG {goal_key} workflow ---")
        results["kegg"] = str(
            run_generic_kegg_target_workflow(
                target_key=goal_key,
                target_config=config,
            )
        )

    elif source_choice == "reactome":
        if "reactome" not in config:
            print(f"{goal_key}에 대한 Reactome workflow 설정이 없습니다.")
            return None

        print(f"\n--- Running Reactome {goal_key} workflow ---")
        reactome_result = run_generic_reactome_target_workflow(
            target_key=goal_key,
            target_config=config,
        )

        if reactome_result is not None:
            results["reactome"] = str(reactome_result)

    elif source_choice == "both":
        if "kegg" in config:
            print(f"\n--- Running KEGG {goal_key} workflow ---")
            results["kegg"] = str(
                run_generic_kegg_target_workflow(
                    target_key=goal_key,
                    target_config=config,
                )
            )

        if "reactome" in config:
            print(f"\n--- Running Reactome {goal_key} workflow ---")
            reactome_result = run_generic_reactome_target_workflow(
                target_key=goal_key,
                target_config=config,
            )

            if reactome_result is not None:
                results["reactome"] = str(reactome_result)

    else:
        return None

    return results or None

def save_combined_workflow_summary(goal_key: str, workflow_results: dict) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_dir = Path("results") / f"{goal_key}_combined_summary_{timestamp}"
    summary_dir.mkdir(parents=True, exist_ok=True)

    summary = {
        "goal": goal_key,
        "workflow_results": workflow_results,
        "created_at": timestamp,
    }

    summary_path = summary_dir / "combined_workflow_summary.json"

    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print(f"\n통합 workflow summary 저장 완료: {summary_path}")
    return summary_path

def run_manual_mode() -> None:
    print(f"--- Manual Mode ---")

    source = input("Enter data source [kegg/reactome]: ").strip().lower()
    pathway = SignalingPathway(source=source)

    try:
        if pathway.source == PathwaySource.KEGG:
            selected_id = _select_kegg_pathway(pathway)

        elif pathway.source == PathwaySource.REACTOME:
            selected_id = _select_reactome_pathway()

        print(f"\n--- Processing {pathway.source.value.upper()} pathway: {selected_id} ---")

        build_reactome_flow = True
        recursive_events = False
        recursive_max_depth = 3

        if pathway.source == PathwaySource.REACTOME:
            flow_answer = input(
                "\nReactome reaction detail(input/output/catalyst)을 확장하시겠습니까? "
                "(y/n, 기본 y): "
            ).strip().lower()

            if flow_answer == "n":
                build_reactome_flow = False

            recursive_answer = input(
                "하위 Pathway의 containedEvents를 재귀적으로 확장하시겠습니까? "
                "(y/n, 기본 n): "
            ).strip().lower()

            if recursive_answer == "y":
                recursive_events = True

                depth_input = input("재귀 확장 depth를 입력하세요. 기본 3: ").strip()
                if depth_input:
                    recursive_max_depth = int(depth_input)

        result = pathway.fetch_pathway(
            selected_id,
            build_reactome_flow=build_reactome_flow,
            recursive_events=recursive_events,
            recursive_max_depth=recursive_max_depth,
        )

        print("\n--- Graph Summary ---")
        print(f"Source: {pathway.graph.graph.get('source')}")
        print(f"Pathway ID: {pathway.graph.graph.get('pathway_id')}")
        print(f"Reactome schemaClass: {pathway.graph.graph.get('reactome_schemaClass', '')}")
        print(f"Mode: {result.get('mode', '')}")
        print(f"Nodes: {len(pathway.graph.nodes)}")
        print(f"Edges: {len(pathway.graph.edges)}")
        print(f"Parser status: {result.get('status')}")

        _print_node_preview(pathway)

        save_summary = input("\nGraph Summary를 별도로 저장하시겠습니까? (y/n): ").strip().lower()
        if save_summary == "y":
            source_name = pathway.graph.graph.get("source", "unknown")
            pathway_id = pathway.graph.graph.get("pathway_id", "unknown")
            safe_pathway_id = str(pathway_id).replace(":", "_").replace("/", "_")

            summary_dir = Path.cwd() / "results" / f"{source_name}_{safe_pathway_id}_summary"
            summary_dir.mkdir(parents=True, exist_ok=True)

            saved_files = pathway.save_graph_summary(
                result=result,
                output_dir=str(summary_dir),
                node_preview_limit=30,
            )

            print("\nGraph Summary 저장 완료:")
            for name, path in saved_files.items():
                print(f"  {name}: {path}")

        show_graph = input("\n그래프를 시각화하시겠습니까? (y/n): ").strip().lower()
        if show_graph == "y":
            source_name = pathway.graph.graph.get("source", "unknown")
            pathway_id = pathway.graph.graph.get("pathway_id", "unknown")
            safe_pathway_id = str(pathway_id).replace(":", "_").replace("/", "_")

            graph_path = Path("results") / f"pathway_graph_{source_name}_{safe_pathway_id}.png"

            pathway.visualize_nodesedges(
                with_edge_labels=False,
                save_path=graph_path,
                show=True,
            )

        _interactive_label_search(pathway)

        paths = _interactive_path_search(pathway)

        _save_all_results(
            pathway=pathway,
            selected_id=selected_id,
            result=result,
            paths=paths,
        )

    except TimeoutError as exc:
        print(f"\n네트워크 시간 초과: {exc}")
        print("잠시 후 다시 시도하거나 timeout 값을 늘려 보세요.")

    except ConnectionError as exc:
        print(f"\n네트워크 오류: {exc}")
        print("인터넷 연결, 방화벽, 프록시 설정을 확인하세요.")

    except Exception as exc:
        import traceback
        print(f"Unexpected Error: {exc}")
        traceback.print_exc()

def select_unregistered_target_mode(target_text: str) -> str:
    print(f"\n'{target_text}'은/는 registry에 등록되어 있지 않습니다.")
    print("진행 방식을 선택하세요.")
    print("  1. Semi-auto registration")
    print("     - KEGG/Reactome 후보 자동 검색")
    print("     - 사용자가 후보 선택")
    print("     - registry에 저장 가능")
    print("  2. Manual mode")
    print("     - KEGG/Reactome ID 직접 입력")
    print("  3. Exit")

    choice = input("선택 [1/2/3]: ").strip().lower()

    if choice in {"1", "semi", "semi-auto", "semi_auto"}:
        return "semi_auto"

    if choice in {"2", "manual"}:
        return "manual"

    return "exit"

def run_semi_auto_registration_flow(target_text: str) -> str | None:
    from workflows.semi_auto_registration import semi_auto_register_target

    try:
        return semi_auto_register_target(target_text)
    except Exception as exc:
        print(f"\nSemi-auto registration 중 오류가 발생했습니다: {exc}")
        return None

def run_post_modeling_analysis(
    target_key: str,
    workflow_results: dict,
    *,
    output_root: str = "results/modeling",
    sbml_output_root: str = "data/sbml/reactome",
    auto_download_sbml: bool = True,
) -> dict:
    """
    KEGG/Reactome workflow 완료 후 모델링 관련 후처리를 수행한다.

    수행 항목:
      1. Reactome nodes.csv 기준 target type classification
      2. Reactome result dir 기준 model evidence check
      3. KEGG nodes/edges 기준 Boolean rule candidate 생성
      4. Reactome SBML 다운로드
      5. SBML checker 실행
      6. SBML summary를 포함한 model evidence 재실행
    """
    print("\n==============================")
    print("Post-modeling analysis 시작")
    print("==============================")

    config = POST_MODELING_CONFIG.get(target_key, {})
    target_keywords = config.get("target_keywords", [])
    reactome_sbml_ids = config.get("reactome_sbml_ids", [])
    sbml_single_check_id = config.get("sbml_single_check_id")

    output_dir = Path(output_root) / target_key
    output_dir.mkdir(parents=True, exist_ok=True)

    post_results: dict = {
        "target_key": target_key,
        "output_dir": str(output_dir),
        "steps": {},
    }

    kegg_dir = workflow_results.get("kegg")
    reactome_dir = workflow_results.get("reactome")

    kegg_dir = Path(kegg_dir) if kegg_dir else None
    reactome_dir = Path(reactome_dir) if reactome_dir else None

    # ------------------------------------------------------------
    # 1. Target type classification
    # ------------------------------------------------------------
    if reactome_dir and (reactome_dir / "nodes.csv").exists() and target_keywords:
        for keyword in target_keywords:
            output_csv = output_dir / f"{keyword}_target_classification.csv"

            cmd = [
                sys.executable,
                "-m",
                _module_name("target_type_classifier"),
                str(reactome_dir / "nodes.csv"),
                "--keyword",
                keyword,
                "--output",
                str(output_csv),
            ]

            ok = _run_command(
                cmd,
                description=f"Target type classification: {keyword}",
            )

            post_results["steps"][f"target_type_{keyword}"] = {
                "status": "success" if ok else "failed",
                "output": str(output_csv),
            }
    else:
        print("[skip] Reactome nodes.csv 또는 target_keywords가 없어 target classification을 건너뜁니다.")

    # ------------------------------------------------------------
    # 2. Model evidence check, SBML 없이 1차 실행
    # ------------------------------------------------------------
    evidence_json = output_dir / "model_evidence.json"
    evidence_md = output_dir / "model_evidence.md"

    if reactome_dir and reactome_dir.exists():
        cmd = [
            sys.executable,
            "-m",
            _module_name("model_evidence_checker"),
            str(reactome_dir),
            "--targets",
            *target_keywords,
            "--output-json",
            str(evidence_json),
            "--output-md",
            str(evidence_md),
        ]

        ok = _run_command(
            cmd,
            description="Model evidence check without SBML",
        )

        post_results["steps"]["model_evidence"] = {
            "status": "success" if ok else "failed",
            "json": str(evidence_json),
            "md": str(evidence_md),
        }
    else:
        print("[skip] Reactome 결과 폴더가 없어 model evidence check를 건너뜁니다.")

    # ------------------------------------------------------------
    # 3. KEGG Boolean model candidate
    # ------------------------------------------------------------
    if (
        kegg_dir
        and (kegg_dir / "nodes.csv").exists()
        and (kegg_dir / "edges.csv").exists()
    ):
        boolean_bnet = output_dir / f"{target_key}_boolean_candidate.bnet"
        boolean_md = output_dir / f"{target_key}_boolean_candidate.md"
        boolean_json = output_dir / f"{target_key}_boolean_candidate.json"

        cmd = [
            sys.executable,
            "-m",
            _module_name("boolean_model_builder"),
            "--nodes",
            str(kegg_dir / "nodes.csv"),
            "--edges",
            str(kegg_dir / "edges.csv"),
            "--include-via-group",
            "--output-json",
            str(boolean_json),
            "--output-bnet",
            str(boolean_bnet),
            "--output-md",
            str(boolean_md),
        ]

        ok = _run_command(
            cmd,
            description="Boolean model candidate from KEGG graph",
        )

        post_results["steps"]["boolean_model_candidate"] = {
            "status": "success" if ok else "failed",
            "json": str(boolean_json),
            "bnet": str(boolean_bnet),
            "md": str(boolean_md),
        }
    else:
        print("[skip] KEGG nodes.csv/edges.csv가 없어 Boolean model candidate 생성을 건너뜁니다.")

    # ------------------------------------------------------------
    # 4. Reactome SBML 다운로드
    # ------------------------------------------------------------
    sbml_target_dir = Path(sbml_output_root) / target_key
    sbml_summary_json = output_dir / f"{target_key}_sbml_directory_summary.json"

    if auto_download_sbml and reactome_sbml_ids:
        cmd = [
            sys.executable,
            "-m",
            _module_name("download_reactome_sbml"),
            "--target-key",
            target_key,
            "--ids",
            *reactome_sbml_ids,
            "--output-root",
            sbml_output_root,
        ]

        ok = _run_command(
            cmd,
            description="Download Reactome SBML",
        )

        post_results["steps"]["download_reactome_sbml"] = {
            "status": "success" if ok else "failed",
            "output_dir": str(sbml_target_dir),
        }
    else:
        print("[skip] Reactome SBML 다운로드 설정이 없어 다운로드를 건너뜁니다.")

    # ------------------------------------------------------------
    # 5. SBML directory checker
    # ------------------------------------------------------------
    if sbml_target_dir.exists():
        cmd = [
            sys.executable,
            "-m",
            _module_name("sbml_checker"),
            str(sbml_target_dir),
            "--output-json",
            str(sbml_summary_json),
        ]

        ok = _run_command(
            cmd,
            description="SBML directory checker",
        )

        post_results["steps"]["sbml_directory_checker"] = {
            "status": "success" if ok else "failed",
            "json": str(sbml_summary_json),
        }
    else:
        print(f"[skip] SBML 디렉토리가 없습니다: {sbml_target_dir}")

    # ------------------------------------------------------------
    # 6. 특정 SBML 파일 1개도 별도 검사
    # ------------------------------------------------------------
    if sbml_single_check_id:
        single_sbml = sbml_target_dir / f"{sbml_single_check_id}.sbml"
        single_summary_json = output_dir / f"{target_key}_{sbml_single_check_id}_sbml_summary.json"

        if single_sbml.exists():
            cmd = [
                sys.executable,
                "-m",
                _module_name("sbml_checker"),
                str(single_sbml),
                "--output-json",
                str(single_summary_json),
            ]

            ok = _run_command(
                cmd,
                description=f"Single SBML checker: {sbml_single_check_id}",
            )

            post_results["steps"]["single_sbml_checker"] = {
                "status": "success" if ok else "failed",
                "json": str(single_summary_json),
            }
        else:
            print(f"[skip] 단일 SBML 파일이 없습니다: {single_sbml}")

    # ------------------------------------------------------------
    # 7. SBML summary를 포함해 model evidence 재실행
    # ------------------------------------------------------------
    evidence_with_sbml_json = output_dir / "model_evidence_with_sbml.json"
    evidence_with_sbml_md = output_dir / "model_evidence_with_sbml.md"

    if reactome_dir and reactome_dir.exists() and sbml_summary_json.exists():
        cmd = [
            sys.executable,
            "-m",
            _module_name("model_evidence_checker"),
            str(reactome_dir),
            "--targets",
            *target_keywords,
            "--sbml-summary-json",
            str(sbml_summary_json),
            "--output-json",
            str(evidence_with_sbml_json),
            "--output-md",
            str(evidence_with_sbml_md),
        ]

        ok = _run_command(
            cmd,
            description="Model evidence check with SBML summary",
        )

        post_results["steps"]["model_evidence_with_sbml"] = {
            "status": "success" if ok else "failed",
            "json": str(evidence_with_sbml_json),
            "md": str(evidence_with_sbml_md),
        }

    # ------------------------------------------------------------
    # 8. post summary 저장
    # ------------------------------------------------------------
    post_summary_path = output_dir / "post_modeling_summary.json"
    with post_summary_path.open("w", encoding="utf-8") as f:
        json.dump(post_results, f, ensure_ascii=False, indent=2)

    print("\nPost-modeling analysis 완료.")
    print(f"Summary: {post_summary_path}")

    return post_results

def main() -> None:
    print(f"--- Environment Check: lxml version {lxml.__version__} ---")

    print("\n분석하고 싶은 Target을 입력하세요.")
    print("예: HIF-1, HIF1A, hypoxia, p53, Keratin, MAPK, PI3K-Akt, VEGF")
    target_text = input("Analysis target: ").strip()

    target_key = recognize_analysis_goal(target_text)

    if target_key == "manual":
        print("\n등록된 자동 workflow가 없습니다.")
        next_mode = select_unregistered_target_mode(target_text)

        if next_mode == "semi_auto":
            target_key = run_semi_auto_registration_flow(target_text)

            if target_key is None:
                print("\nSemi-auto registration이 완료되지 않았습니다. Manual mode로 전환합니다.")
                run_manual_mode()
                return

        elif next_mode == "manual":
            run_manual_mode()
            return

        else:
            print("종료합니다.")
            return

    print_goal_workflow_description(target_key)

    source_choice = select_data_source_for_goal()

    if source_choice == "manual":
        run_manual_mode()
        return

    workflow_results = run_goal_workflow(target_key, source_choice)

    if workflow_results:
        save_combined_workflow_summary(target_key, workflow_results)

        print("\nWorkflow 완료.")
        for source_name, output_dir in workflow_results.items():
            print(f"  {source_name}: {output_dir}")

        run_post = input("\n모델링 후처리 분석을 실행하시겠습니까? (y/n): ").strip().lower()

        if run_post == "y":
            run_post_modeling_analysis(
                target_key=target_key,
                workflow_results=workflow_results,
                output_root="results/modeling",
                sbml_output_root="data/sbml/reactome",
                auto_download_sbml=True,
            )

    else:
        print("\nworkflow 실행 결과가 없습니다. Manual mode로 전환합니다.")
        run_manual_mode()
  
if __name__ == "__main__":
    main()
    
