# -*- coding: utf-8 -*-
"""
Created on Mon May 25 13:17:47 2026

@author: user
"""

from pathlib import Path
import pandas as pd
import networkx as nx

result_dir = Path(r"results\kegg_hif1_hsa04066_YYYYMMDD_HHMMSS")

nodes = pd.read_csv(result_dir / "nodes.csv")
edges = pd.read_csv(result_dir / "edges.csv")

graph = nx.DiGraph()

for _, row in nodes.iterrows():
    graph.add_node(
        row["node_id"],
        label=row.get("label", ""),
        type=row.get("type", ""),
    )

for _, row in edges.iterrows():
    graph.add_edge(
        row["source_node"],
        row["target_node"],
        relation_type=row.get("relation_type", ""),
    )

start = "KEGG:3091"
end = "KEGG:7422"

print("start exists:", start in graph)
print("end exists:", end in graph)

if start in graph:
    print("start label:", graph.nodes[start].get("label"))

if end in graph:
    print("end label:", graph.nodes[end].get("label"))

print("directed path HIF1A -> VEGFA:", nx.has_path(graph, start, end) if start in graph and end in graph else False)
print("reverse path VEGFA -> HIF1A:", nx.has_path(graph, end, start) if start in graph and end in graph else False)
print("undirected connection:", nx.has_path(graph.to_undirected(), start, end) if start in graph and end in graph else False)
