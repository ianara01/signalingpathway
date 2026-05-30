# -*- coding: utf-8 -*-
"""Parser for BlastKOALA / GhostKOALA / KofamKOALA-like outputs."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import re
import pandas as pd

KO_PATTERN = re.compile(r"\bK\d{5}\b")


def parse_koala_table(text: str, source_format: str = "auto") -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        significant = False
        if line.startswith("*"):
            significant = True
            line = line[1:].strip()
        sep = "\t" if "\t" in line else ("," if "," in line else None)
        parts = [p.strip() for p in line.split(sep)] if sep else line.split()
        match = KO_PATTERN.search(line)
        if not match:
            continue
        ko_id = match.group(0)
        query_id = parts[0] if parts else ""
        annotation = ""
        if len(parts) >= 3:
            ko_index = next((i for i, p in enumerate(parts) if KO_PATTERN.fullmatch(p)), None)
            if ko_index is not None:
                annotation = " ".join(parts[ko_index + 1:]).strip()
        numeric = [p for p in parts if re.fullmatch(r"[-+]?\d*\.?\d+(e[-+]?\d+)?", p, re.I)]
        rows.append({
            "query_id": query_id,
            "ko_id": ko_id,
            "score": numeric[-1] if numeric else "",
            "threshold": numeric[0] if len(numeric) >= 2 else "",
            "e_value": next((p for p in parts if re.fullmatch(r"[-+]?\d*\.?\d+e[-+]?\d+", p, re.I)), ""),
            "significant": significant,
            "annotation": annotation,
            "source_format": source_format,
            "raw_line": raw,
        })
    return pd.DataFrame(rows)


def parse_koala_file(path: str | Path, source_format: str = "auto") -> pd.DataFrame:
    return parse_koala_table(Path(path).read_text(encoding="utf-8", errors="ignore"), source_format=source_format)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Parse KOALA output.")
    parser.add_argument("input_file")
    parser.add_argument("--output", default=None)
    parser.add_argument("--source-format", default="auto")
    args = parser.parse_args()
    df = parse_koala_file(args.input_file, args.source_format)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, index=False, encoding="utf-8-sig")
        print(f"Saved: {args.output}")
    else:
        print(df.to_string(index=False))


if __name__ == "__main__":
    main()
