# -*- coding: utf-8 -*-
"""Extended KEGG provider.

Adds KEGG GENES / COMPOUND / DRUG / DISEASE / ENZYME / REACTION coverage helpers
using KEGG REST operations: get, list, find, conv, link.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import time
import requests
import pandas as pd


class KeggExtendedProvider:
    BASE_URL = "https://rest.kegg.jp"

    def __init__(self, timeout: int = 30, sleep_sec: float = 0.2):
        self.timeout = timeout
        self.sleep_sec = sleep_sec
        self._cache: dict[str, str] = {}

    def _request_text(self, url: str) -> str:
        if url in self._cache:
            return self._cache[url]
        response = requests.get(url, timeout=self.timeout, headers={"User-Agent": "S2SignalPathway/1.0"})
        response.raise_for_status()
        text = response.text
        self._cache[url] = text
        time.sleep(self.sleep_sec)
        return text

    def get(self, entry: str) -> str:
        return self._request_text(f"{self.BASE_URL}/get/{entry}")

    def list_db(self, database: str) -> str:
        return self._request_text(f"{self.BASE_URL}/list/{database}")

    def find(self, database: str, query: str) -> str:
        return self._request_text(f"{self.BASE_URL}/find/{database}/{query}")

    def link(self, target_db: str, source_db_or_entry: str) -> str:
        return self._request_text(f"{self.BASE_URL}/link/{target_db}/{source_db_or_entry}")

    def conv(self, target_db: str, source_db_or_entry: str) -> str:
        return self._request_text(f"{self.BASE_URL}/conv/{target_db}/{source_db_or_entry}")

    @staticmethod
    def parse_kegg_flat_entry(text: str) -> dict[str, list[str]]:
        fields: dict[str, list[str]] = {}
        current_key = None
        for line in text.splitlines():
            if not line.strip():
                continue
            key = line[:12].strip()
            value = line[12:].rstrip()
            if key:
                current_key = key
                fields.setdefault(current_key, []).append(value.strip())
            elif current_key:
                fields[current_key].append(value.strip())
        return fields

    @staticmethod
    def first_field(fields: dict[str, list[str]], key: str) -> str:
        values = fields.get(key, [])
        return values[0] if values else ""

    @staticmethod
    def join_field(fields: dict[str, list[str]], key: str, sep: str = "; ") -> str:
        return sep.join([x for x in fields.get(key, []) if x])

    @staticmethod
    def parse_link_table(text: str) -> pd.DataFrame:
        rows = []
        for line in text.splitlines():
            parts = line.strip().split("\t")
            if len(parts) == 2:
                rows.append({"source": parts[0], "target": parts[1]})
        return pd.DataFrame(rows)

    @staticmethod
    def normalize_hsa_gene_id(value: str) -> str:
        v = str(value).strip().replace("KEGG:", "")
        if v.startswith("hsa:"):
            return v
        if v.isdigit():
            return f"hsa:{v}"
        if v.startswith("hsa") and ":" not in v:
            return v.replace("hsa", "hsa:", 1)
        return v

    @staticmethod
    def normalize_compound_id(value: str) -> str:
        return str(value).strip().replace("KEGG:", "").replace("cpd:", "").replace("compound:", "")

    def annotate_gene(self, gene_id: str) -> dict[str, Any]:
        entry = self.normalize_hsa_gene_id(gene_id)
        fields = self.parse_kegg_flat_entry(self.get(entry))
        name = self.first_field(fields, "NAME")
        return {
            "kegg_gene_id": entry,
            "entry": self.first_field(fields, "ENTRY"),
            "name": name,
            "symbol": name.split(",")[0].strip() if name else "",
            "definition": self.first_field(fields, "DEFINITION"),
            "organism": self.first_field(fields, "ORGANISM"),
            "orthology": self.join_field(fields, "ORTHOLOGY"),
            "pathway": self.join_field(fields, "PATHWAY"),
            "brite": self.join_field(fields, "BRITE"),
            "disease": self.join_field(fields, "DISEASE"),
            "drug_target": self.join_field(fields, "DRUG_TARGET"),
            "dblinks": self.join_field(fields, "DBLINKS"),
            "aa_seq_present": "AASEQ" in fields,
            "nt_seq_present": "NTSEQ" in fields,
        }

    def annotate_compound(self, compound_id: str) -> dict[str, Any]:
        entry = self.normalize_compound_id(compound_id)
        fields = self.parse_kegg_flat_entry(self.get(entry))
        return {
            "compound_id": entry,
            "entry": self.first_field(fields, "ENTRY"),
            "name": self.join_field(fields, "NAME"),
            "formula": self.first_field(fields, "FORMULA"),
            "exact_mass": self.first_field(fields, "EXACT_MASS"),
            "mol_weight": self.first_field(fields, "MOL_WEIGHT"),
            "reaction": self.join_field(fields, "REACTION"),
            "pathway": self.join_field(fields, "PATHWAY"),
            "enzyme": self.join_field(fields, "ENZYME"),
            "brite": self.join_field(fields, "BRITE"),
            "dblinks": self.join_field(fields, "DBLINKS"),
        }

    def annotate_drug(self, drug_id: str) -> dict[str, Any]:
        fields = self.parse_kegg_flat_entry(self.get(drug_id))
        return {
            "drug_id": drug_id,
            "entry": self.first_field(fields, "ENTRY"),
            "name": self.join_field(fields, "NAME"),
            "formula": self.first_field(fields, "FORMULA"),
            "target": self.join_field(fields, "TARGET"),
            "pathway": self.join_field(fields, "PATHWAY"),
            "disease": self.join_field(fields, "DISEASE"),
            "remark": self.join_field(fields, "REMARK"),
            "dblinks": self.join_field(fields, "DBLINKS"),
        }

    def annotate_disease(self, disease_id: str) -> dict[str, Any]:
        fields = self.parse_kegg_flat_entry(self.get(disease_id))
        return {
            "disease_id": disease_id,
            "entry": self.first_field(fields, "ENTRY"),
            "name": self.join_field(fields, "NAME"),
            "description": self.join_field(fields, "DESCRIPTION"),
            "category": self.join_field(fields, "CATEGORY"),
            "gene": self.join_field(fields, "GENE"),
            "pathway": self.join_field(fields, "PATHWAY"),
            "drug": self.join_field(fields, "DRUG"),
            "dblinks": self.join_field(fields, "DBLINKS"),
        }

    def annotate_enzyme(self, enzyme_id: str) -> dict[str, Any]:
        fields = self.parse_kegg_flat_entry(self.get(enzyme_id))
        return {
            "enzyme_id": enzyme_id,
            "entry": self.first_field(fields, "ENTRY"),
            "name": self.join_field(fields, "NAME"),
            "class": self.join_field(fields, "CLASS"),
            "sysname": self.join_field(fields, "SYSNAME"),
            "reaction": self.join_field(fields, "REACTION"),
            "substrate": self.join_field(fields, "SUBSTRATE"),
            "product": self.join_field(fields, "PRODUCT"),
            "pathway": self.join_field(fields, "PATHWAY"),
            "genes": self.join_field(fields, "GENES"),
            "dblinks": self.join_field(fields, "DBLINKS"),
        }

    def annotate_reaction(self, reaction_id: str) -> dict[str, Any]:
        fields = self.parse_kegg_flat_entry(self.get(reaction_id))
        return {
            "reaction_id": reaction_id,
            "entry": self.first_field(fields, "ENTRY"),
            "name": self.join_field(fields, "NAME"),
            "definition": self.join_field(fields, "DEFINITION"),
            "equation": self.join_field(fields, "EQUATION"),
            "enzyme": self.join_field(fields, "ENZYME"),
            "pathway": self.join_field(fields, "PATHWAY"),
            "rclass": self.join_field(fields, "RCLASS"),
        }

    def link_entry_to_databases(self, entry: str, target_dbs: list[str]) -> dict[str, pd.DataFrame]:
        out = {}
        for db in target_dbs:
            try:
                out[db] = self.parse_link_table(self.link(db, entry))
            except Exception as exc:
                out[db] = pd.DataFrame([{"source": entry, "target": "", "error": str(exc)}])
        return out

    def build_gene_coverage(self, gene_ids: list[str]) -> dict[str, pd.DataFrame]:
        records, link_tables = [], []
        for gid in gene_ids:
            entry = self.normalize_hsa_gene_id(gid)
            try:
                records.append(self.annotate_gene(entry))
            except Exception as exc:
                records.append({"kegg_gene_id": entry, "error": str(exc)})
            for db, df in self.link_entry_to_databases(entry, ["pathway", "disease", "drug", "enzyme", "reaction", "ko"]).items():
                if not df.empty:
                    df = df.copy()
                    df["target_database"] = db
                    link_tables.append(df)
        return {
            "gene_annotations": pd.DataFrame(records),
            "links": pd.concat(link_tables, ignore_index=True) if link_tables else pd.DataFrame(),
        }

    def build_compound_coverage(self, compound_ids: list[str]) -> dict[str, pd.DataFrame]:
        records, link_tables = [], []
        for cid in compound_ids:
            entry = self.normalize_compound_id(cid)
            try:
                records.append(self.annotate_compound(entry))
            except Exception as exc:
                records.append({"compound_id": entry, "error": str(exc)})
            for db, df in self.link_entry_to_databases(entry, ["pathway", "reaction", "enzyme", "drug"]).items():
                if not df.empty:
                    df = df.copy()
                    df["target_database"] = db
                    link_tables.append(df)
        return {
            "compound_annotations": pd.DataFrame(records),
            "links": pd.concat(link_tables, ignore_index=True) if link_tables else pd.DataFrame(),
        }


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Extended KEGG coverage.")
    parser.add_argument("--genes", nargs="*", default=[])
    parser.add_argument("--compounds", nargs="*", default=[])
    parser.add_argument("--output-dir", default="results/kegg_extended")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    provider = KeggExtendedProvider()

    summary = {"genes": args.genes, "compounds": args.compounds, "output_dir": str(out)}
    if args.genes:
        cov = provider.build_gene_coverage(args.genes)
        cov["gene_annotations"].to_csv(out / "kegg_gene_annotations.csv", index=False, encoding="utf-8-sig")
        cov["links"].to_csv(out / "kegg_gene_database_links.csv", index=False, encoding="utf-8-sig")
    if args.compounds:
        cov = provider.build_compound_coverage(args.compounds)
        cov["compound_annotations"].to_csv(out / "kegg_compound_annotations.csv", index=False, encoding="utf-8-sig")
        cov["links"].to_csv(out / "kegg_compound_database_links.csv", index=False, encoding="utf-8-sig")
    (out / "kegg_extended_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
