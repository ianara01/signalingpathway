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

# S2SignalPathway extended modules

추가 파일:
- providers/kegg_extended_provider.py
- providers/reactome_analysis_provider.py
- providers/reactome_species_provider.py
- annotations/koala_parser.py
- annotations/ko_to_pathway_mapper.py
- analysis/kegg_database_coverage.py
- analysis/reactome_enrichment.py
- analysis/expression_mapping.py
- analysis/species_comparison.py
- main_patch/main_extended_analysis_patch.py

## main.py 연결

`main_patch/main_extended_analysis_patch.py`의 import/function들을 `main.py`에 붙인 뒤,
기존 `main()`에서 workflow 완료 후 아래 블록을 추가합니다.

```python
run_ext = input("\nKEGG/Reactome 확장 분석을 실행하시겠습니까? (y/n): ").strip().lower()
if run_ext == "y":
    koala_file = input("KOALA 결과 파일 경로, 없으면 Enter: ").strip() or None
    expression_file = input("Expression CSV 파일 경로, 없으면 Enter: ").strip() or None
    run_extended_database_analysis(
        target_key=target_key,
        workflow_results=workflow_results,
        output_root="results/analysis_extended",
        koala_file=koala_file,
        expression_file=expression_file,
    )
```

## 독립 실행

```bat
python -m analysis.kegg_database_coverage results\kegg_hif1_hsa04066_YYYYMMDD_HHMMSS --output-dir results\analysis_extended\hif1\kegg_database_coverage

python -m analysis.reactome_enrichment --identifiers HIF1A VEGFA ARNT VHL EGLN1 EGLN2 --output-dir results\analysis_extended\hif1\reactome_enrichment

python -m annotations.koala_parser data\annotations\ghostkoala.txt --output results\analysis_extended\hif1\annotations\koala_parsed.csv

python -m annotations.ko_to_pathway_mapper results\analysis_extended\hif1\annotations\koala_parsed.csv --output-csv results\analysis_extended\hif1\annotations\ko_pathway_mapping.csv

python -m analysis.expression_mapping --nodes results\reactome_hif1_combined_YYYYMMDD_HHMMSS\nodes.csv --expression data\expression\expression_matrix.csv --output-dir results\analysis_extended\hif1\expression_mapping

python -m analysis.species_comparison --reactome-result-dir results\reactome_hif1_combined_YYYYMMDD_HHMMSS --species "Mus musculus" "Rattus norvegicus" --output-dir results\analysis_extended\hif1\species_comparison
```

## 주요 출력

```text
results/analysis_extended/hif1/
├─ kegg_database_coverage/
│  ├─ kegg_database_coverage_matrix.csv
│  ├─ kegg_gene_annotations.csv
│  ├─ kegg_gene_database_links.csv
│  ├─ kegg_compound_annotations.csv
│  └─ kegg_database_coverage_summary.json
├─ reactome_enrichment/
├─ annotations/
├─ expression_mapping/
├─ species_comparison/
└─ extended_analysis_summary.json
```
