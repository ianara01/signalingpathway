# -*- coding: utf-8 -*-
"""
Created on Thu May 21 00:26:23 2026

@author: user
"""

from pathlib import Path
import re
import pandas as pd
import networkx as nx


def sanitize_sheet_name(name: str, max_len: int = 31) -> str:
    """
    Excel sheet name 제한 처리:
    - 최대 31자
    - 특수문자 []:*?/\\ 제거
    """
    name = re.sub(r"[\[\]\:\*\?\/\\]", "_", str(name))
    return name[:max_len]


def load_nodes(nodes_csv_path: Path) -> dict:
    """
    nodes.csv를 읽어서 node_id별 속성 dict 생성.
    """
    if not nodes_csv_path.exists():
        raise FileNotFoundError(f"nodes.csv not found: {nodes_csv_path}")

    nodes_df = pd.read_csv(nodes_csv_path)

    required = {"node_id"}
    missing = required - set(nodes_df.columns)
    if missing:
        raise ValueError(f"nodes.csv에 필수 컬럼이 없습니다: {missing}")

    node_attrs = {}

    for _, row in nodes_df.iterrows():
        node_id = str(row.get("node_id", "")).strip()
        if not node_id:
            continue

        node_attrs[node_id] = {
            "source": row.get("source", ""),
            "type": row.get("type", ""),
            "label": row.get("label", ""),
            "raw_id": row.get("raw_id", ""),
        }

    return node_attrs


def build_graph_from_nodes_edges(
    nodes_csv_path: Path,
    edges_csv_path: Path,
) -> tuple[nx.DiGraph, pd.DataFrame]:
    """
    nodes.csv와 edges.csv를 읽어 directed graph 생성.
    """
    if not edges_csv_path.exists():
        raise FileNotFoundError(f"edges.csv not found: {edges_csv_path}")

    node_attrs = load_nodes(nodes_csv_path)
    edges_df = pd.read_csv(edges_csv_path)

    required = {"source_node", "target_node"}
    missing = required - set(edges_df.columns)
    if missing:
        raise ValueError(f"edges.csv에 필수 컬럼이 없습니다: {missing}")

    graph = nx.DiGraph()

    # nodes.csv 기반 노드 추가
    for node_id, attrs in node_attrs.items():
        graph.add_node(
            node_id,
            source=attrs.get("source", ""),
            type=attrs.get("type", ""),
            label=attrs.get("label", ""),
            raw_id=attrs.get("raw_id", ""),
        )

    # edges.csv 기반 edge 추가
    for edge_index, row in edges_df.iterrows():
        source_node = str(row["source_node"]).strip()
        target_node = str(row["target_node"]).strip()

        if source_node not in graph:
            graph.add_node(
                source_node,
                source=row.get("source", ""),
                type="",
                label="",
                raw_id="",
            )

        if target_node not in graph:
            graph.add_node(
                target_node,
                source=row.get("source", ""),
                type="",
                label="",
                raw_id="",
            )

        graph.add_edge(
            source_node,
            target_node,
            relation_type=row.get("relation_type", ""),
            source=row.get("source", ""),
            original_edge_index=edge_index + 1,
        )

    return graph, edges_df


def paths_to_long_dataframe(
    graph: nx.DiGraph,
    paths: list[list[str]],
    edge_pair_index: int,
    start_node: str,
    end_node: str,
    base_relation_type: str,
) -> pd.DataFrame:
    """
    여러 path를 long-format으로 저장.
    paths_index.csv를 세로로 확장한 형태.
    """
    rows = []

    for path_index, path in enumerate(paths, start=1):
        for step_index, node_id in enumerate(path, start=1):
            attrs = graph.nodes[node_id]

            rows.append(
                {
                    "edge_pair_index": edge_pair_index,
                    "path_index": path_index,
                    "step_index": step_index,
                    "start_node": start_node,
                    "end_node": end_node,
                    "base_relation_type": base_relation_type,
                    "node_id": node_id,
                    "label": attrs.get("label", ""),
                    "type": attrs.get("type", ""),
                    "source": attrs.get("source", ""),
                    "raw_id": attrs.get("raw_id", ""),
                }
            )

    return pd.DataFrame(
        rows,
        columns=[
            "edge_pair_index",
            "path_index",
            "step_index",
            "start_node",
            "end_node",
            "base_relation_type",
            "node_id",
            "label",
            "type",
            "source",
            "raw_id",
        ],
    )


def paths_to_index_style_dataframe(
    graph: nx.DiGraph,
    paths: list[list[str]],
) -> pd.DataFrame:
    """
    paths_index.csv와 유사한 가로형 구조.
    하나의 sheet 안에 path별 block을 이어붙인다.

    예:
      path_index | 1
      step_index | 1 | 2 | 3
      node_id    | A | B | C
      label      | ...
      type       | ...
    """
    blocks = []

    for path_index, path in enumerate(paths, start=1):
        max_step = len(path)
        columns = ["field"] + [str(i) for i in range(1, max_step + 1)]

        block_rows = []

        block_rows.append(["path_index"] + [path_index] + [""] * (max_step - 1))
        block_rows.append(["step_index"] + list(range(1, max_step + 1)))
        block_rows.append(["node_id"] + path)
        block_rows.append(
            ["label"] + [graph.nodes[node].get("label", "") for node in path]
        )
        block_rows.append(
            ["type"] + [graph.nodes[node].get("type", "") for node in path]
        )
        block_rows.append(
            ["source"] + [graph.nodes[node].get("source", "") for node in path]
        )
        block_rows.append(
            ["raw_id"] + [graph.nodes[node].get("raw_id", "") for node in path]
        )

        block_df = pd.DataFrame(block_rows, columns=columns)

        # path block 사이 빈 줄
        empty_df = pd.DataFrame([[""] * len(columns)], columns=columns)

        blocks.append(block_df)
        blocks.append(empty_df)

    if not blocks:
        return pd.DataFrame(
            [{"message": "No directed path found"}]
        )

    return pd.concat(blocks, ignore_index=True)


def path_edges_dataframe(
    graph: nx.DiGraph,
    paths: list[list[str]],
    edge_pair_index: int,
) -> pd.DataFrame:
    """
    각 path의 edge 상세 정보 저장.
    """
    rows = []

    for path_index, path in enumerate(paths, start=1):
        for step_index, (u, v) in enumerate(zip(path[:-1], path[1:]), start=1):
            edge_attrs = graph.get_edge_data(u, v, default={})

            rows.append(
                {
                    "edge_pair_index": edge_pair_index,
                    "path_index": path_index,
                    "edge_step_index": step_index,
                    "source_node": u,
                    "source_label": graph.nodes[u].get("label", ""),
                    "target_node": v,
                    "target_label": graph.nodes[v].get("label", ""),
                    "relation_type": edge_attrs.get("relation_type", ""),
                    "source": edge_attrs.get("source", ""),
                }
            )

    return pd.DataFrame(
        rows,
        columns=[
            "edge_pair_index",
            "path_index",
            "edge_step_index",
            "source_node",
            "source_label",
            "target_node",
            "target_label",
            "relation_type",
            "source",
        ],
    )


def find_paths_for_all_edges_to_excel(
    result_dir: str,
    cutoff: int = 10,
    max_paths_per_edge: int | None = None,
    output_filename: str = "paths_index_all_edges.xlsx",
) -> Path:
    """
    같은 디렉토리의 nodes.csv, edges.csv를 읽어서
    edges.csv의 모든 source_node -> target_node pair에 대해
    가능한 모든 directed path를 찾고,
    edge pair별 별도 sheet에 저장한다.
    """
    result_dir = Path(result_dir)

    nodes_csv_path = result_dir / "nodes.csv"
    edges_csv_path = result_dir / "edges.csv"

    graph, edges_df = build_graph_from_nodes_edges(
        nodes_csv_path=nodes_csv_path,
        edges_csv_path=edges_csv_path,
    )

    output_path = result_dir / output_filename

    summary_rows = []
    all_long_rows = []
    all_edge_rows = []

    print(f"Result dir : {result_dir}")
    print(f"nodes.csv  : {nodes_csv_path}")
    print(f"edges.csv  : {edges_csv_path}")
    print(f"Graph nodes: {graph.number_of_nodes()}")
    print(f"Graph edges: {graph.number_of_edges()}")
    print(f"Edge pairs : {len(edges_df)}")

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for edge_pair_index, row in edges_df.iterrows():
            pair_no = edge_pair_index + 1

            start_node = str(row["source_node"]).strip()
            end_node = str(row["target_node"]).strip()
            relation_type = row.get("relation_type", "")
            source = row.get("source", "")

            if start_node not in graph or end_node not in graph:
                paths = []
            elif nx.has_path(graph, start_node, end_node):
                path_generator = nx.all_simple_paths(
                    graph,
                    source=start_node,
                    target=end_node,
                    cutoff=cutoff,
                )

                if max_paths_per_edge is None:
                    paths = list(path_generator)
                else:
                    paths = []
                    for path in path_generator:
                        paths.append(path)
                        if len(paths) >= max_paths_per_edge:
                            break
            else:
                paths = []

            start_label = graph.nodes[start_node].get("label", "") if start_node in graph else ""
            end_label = graph.nodes[end_node].get("label", "") if end_node in graph else ""

            summary_rows.append(
                {
                    "edge_pair_index": pair_no,
                    "start_node": start_node,
                    "start_label": start_label,
                    "end_node": end_node,
                    "end_label": end_label,
                    "base_relation_type": relation_type,
                    "source": source,
                    "path_count": len(paths),
                    "cutoff": cutoff,
                }
            )

            # edge pair별 sheet 저장
            short_start = start_node.replace("KEGG:", "").replace("REACTOME:", "")
            short_end = end_node.replace("KEGG:", "").replace("REACTOME:", "")

            sheet_name = sanitize_sheet_name(
                f"E{pair_no:03d}_{short_start}_to_{short_end}"
            )

            index_style_df = paths_to_index_style_dataframe(graph, paths)
            index_style_df.to_excel(writer, sheet_name=sheet_name, index=False)

            # 전체 long-format 누적
            long_df = paths_to_long_dataframe(
                graph=graph,
                paths=paths,
                edge_pair_index=pair_no,
                start_node=start_node,
                end_node=end_node,
                base_relation_type=relation_type,
            )

            if not long_df.empty:
                all_long_rows.append(long_df)

            edge_path_df = path_edges_dataframe(
                graph=graph,
                paths=paths,
                edge_pair_index=pair_no,
            )

            if not edge_path_df.empty:
                all_edge_rows.append(edge_path_df)

            print(
                f"[{pair_no:03d}] {start_node} -> {end_node} | "
                f"paths={len(paths)} | relation={relation_type}"
            )

        # Summary sheet
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_excel(writer, sheet_name="Summary", index=False)

        # 모든 path를 하나의 long-format sheet에도 저장
        if all_long_rows:
            all_paths_long_df = pd.concat(all_long_rows, ignore_index=True)
        else:
            all_paths_long_df = pd.DataFrame(
                [{"message": "No directed paths found"}]
            )

        all_paths_long_df.to_excel(writer, sheet_name="All_Paths_Long", index=False)

        # 모든 path edge 상세
        if all_edge_rows:
            all_path_edges_df = pd.concat(all_edge_rows, ignore_index=True)
        else:
            all_path_edges_df = pd.DataFrame(
                [{"message": "No path edges found"}]
            )

        all_path_edges_df.to_excel(writer, sheet_name="All_Path_Edges", index=False)

        # column width 조정
        for sheet_name, worksheet in writer.sheets.items():
            worksheet.freeze_panes = "A2"
            for col in range(1, 15):
                worksheet.column_dimensions[chr(64 + col)].width = 22

    print(f"\nSaved Excel workbook: {output_path}")
    return output_path


if __name__ == "__main__":
    result_dir = r"C:\workspaces\S2SignalPathway\results\kegg_hsa04066_20260520_225848"

    find_paths_for_all_edges_to_excel(
        result_dir=result_dir,
        cutoff=10,
        max_paths_per_edge=None,
        output_filename="paths_index_all_edges.xlsx",
    )