# S2SignalPathway Refactored

KEGG와 Reactome 처리를 명확히 분리한 리팩터링 버전입니다.

## 구조

```text
S2SignalPathway/
├─ __init__.py
├─ s2pathway.py
├─ main.py
├─ providers/
│  ├─ __init__.py
│  ├─ kegg_provider.py
│  └─ reactome_provider.py
├─ parsers/
│  ├─ __init__.py
│  ├─ kegg_kgml_parser.py
│  └─ reactome_json_parser.py
└─ visualization/
   ├─ __init__.py
   └─ pathway_visualizer.py
```

## 핵심 변경점

- `SignalingPathway.fetch_pathway()` 단일 진입점 추가
- `source='kegg'`일 때 KEGG ID/API/KGML 파서만 사용
- `source='reactome'`일 때 Reactome ID/API/JSON 파서만 사용
- 노드 ID prefix 적용
  - KEGG: `KEGG:3091`
  - Reactome: `REACTOME:R-HSA-...`
- 엣지 속성명 `relation_type`으로 통일
- 기존 자동 실행되던 `PathwayVisualizer("hsa04010")` 예제 제거

## 실행

```bash
pip install -r requirements.txt
python -m S2SignalPathway.main
```

## 사용 예

```python
from S2SignalPathway import SignalingPathway

# KEGG
kegg = SignalingPathway(source="kegg")
kegg.fetch_pathway("hsa04010")
print(kegg.graph.nodes)

# Reactome
reactome = SignalingPathway(source="reactome")
reactome.fetch_pathway("R-HSA-168256")
print(reactome.graph.nodes)
```
