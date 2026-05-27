"""KEGG API provider.

This module is intentionally KEGG-only. It returns raw KEGG API payloads and
keeps KGML parsing in ``parsers/kegg_kgml_parser.py``.
"""

from __future__ import annotations

import re
from typing import Any

import requests


class KeggProvider:
    """Thin client for KEGG REST API."""

    BASE_URL = "https://rest.kegg.jp"

    def __init__(self, timeout: int = 30):
        self.timeout = timeout
        self._compound_name_cache: dict[str, str | None] = {}

    @staticmethod
    def normalize_pathway_id(pathway_id: str) -> str:
        """Normalize KEGG pathway IDs such as ``path:hsa04010`` or ``HSA04010``."""
        if not pathway_id:
            raise ValueError("KEGG pathway_id is empty.")

        clean_id = pathway_id.strip().replace("path:", "").lower()
        if not clean_id.startswith(("hsa", "mmu", "rno", "map")):
            raise ValueError(
                "KEGG pathway ID 형식이 아닙니다. 예: hsa04010, hsa04151, map04010"
            )
        return clean_id

    def search_pathway_id(self, keyword: str) -> list[dict[str, str]]:
        """Search KEGG pathway IDs by keyword.

        Includes a small convenience search that also tries a hyphenated version
        for terms such as HIF1 -> HIF-1.
        """
        keyword = keyword.strip()
        if not keyword:
            return []

        keywords = [keyword]
        if "-" not in keyword and any(char.isdigit() for char in keyword):
            keywords.append(re.sub(r"(\d+)", r"-\1", keyword))

        all_results: list[dict[str, str]] = []
        for kw in keywords:
            url = f"{self.BASE_URL}/find/pathway/{kw}"
            response = requests.get(url, timeout=self.timeout)
            if response.status_code != 200 or not response.text.strip():
                continue

            for line in response.text.strip().split("\n"):
                parts = line.split("\t", maxsplit=1)
                if len(parts) != 2:
                    continue
                path_id = parts[0].replace("path:", "")
                description = parts[1]
                all_results.append({"id": path_id, "name": description})

        return list({item["id"]: item for item in all_results}.values())

    def fetch_kgml(self, pathway_id: str) -> str:
        """Fetch KEGG KGML XML for one pathway."""
        clean_id = self.normalize_pathway_id(pathway_id)
        url = f"{self.BASE_URL}/get/{clean_id}/kgml"
        response = requests.get(url, timeout=self.timeout)
        if response.status_code != 200 or not response.text.strip():
            raise ValueError(f"Failed to fetch KEGG KGML data for {clean_id}")
        return response.text

    def fetch_pathway_payload(self, pathway_id: str) -> dict[str, Any]:
        """Return a provider-neutral payload consumed by the KEGG parser."""
        clean_id = self.normalize_pathway_id(pathway_id)
        return {
            "source": "kegg",
            "pathway_id": clean_id,
            "format": "kgml",
            "data": self.fetch_kgml(clean_id),
        }

    def normalize_compound_id(self, compound_id: str) -> str:
        """
        KEGG compound ID를 Cxxxxx 형식으로 정규화한다.
    
        허용 예:
          C00007
          KEGG:C00007
          cpd:C00007
          compound:C00007
        """
        clean_id = (
            str(compound_id)
            .strip()
            .replace("KEGG:", "")
            .replace("kegg:", "")
            .replace("cpd:", "")
            .replace("CPD:", "")
            .replace("compound:", "")
            .replace("COMPOUND:", "")
        )
    
        if not clean_id.startswith("C"):
            raise ValueError(f"Invalid KEGG compound ID: {compound_id}")
    
        return clean_id
    
    def fetch_compound_entry(self, compound_id: str) -> str:
        clean_id = self.normalize_compound_id(compound_id)
    
        url = f"{self.BASE_URL}/get/{clean_id}"
        response = requests.get(url, timeout=self.timeout)
        response.raise_for_status()
        return response.text

    def fetch_compound_name(self, compound_id: str) -> str | None:
        """
        KEGG compound ID에서 NAME 필드를 추출한다.
        """
        try:
            clean_id = self.normalize_compound_id(compound_id)
        except ValueError:
            return None
    
        if clean_id in self._compound_name_cache:
            return self._compound_name_cache[clean_id]
    
        try:
            text = self.fetch_compound_entry(clean_id)
        except Exception:
            self._compound_name_cache[clean_id] = None
            return None
    
        names: list[str] = []
        in_name = False
    
        for line in text.splitlines():
            if line.startswith("NAME"):
                in_name = True
                value = line[12:].strip()
                if value:
                    names.append(value.rstrip(";"))
    
            elif in_name and line.startswith(" " * 12):
                value = line[12:].strip()
                if value:
                    names.append(value.rstrip(";"))
    
            elif in_name:
                break
    
        if not names:
            self._compound_name_cache[clean_id] = None
            return None
    
        self._compound_name_cache[clean_id] = names[0]
        return names[0]
