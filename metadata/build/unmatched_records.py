#!/usr/bin/env python3
"""List maize cold RNA-seq repository records that match no catalogue study (for Rubén to triage).

Inputs: metadata/build/repo_sra_maize_cold.tsv and repo_geo_maize_cold.tsv (search_repositories.py),
metadata/accessions_resolved.tsv and metadata/build/accession_records.tsv (linked IDs).
Output: metadata/build/uncatalogued_repository_records.tsv. Nothing is added to the catalogue.
Usage (project root): python3 metadata/build/unmatched_records.py
"""
import csv
import re
from pathlib import Path

BUILD = Path(__file__).resolve().parent
known = set()
for r in csv.DictReader(open(BUILD.parent / "accessions_resolved.tsv"), delimiter="\t"):
    known.add(r["accession"])
    known.update(re.findall(r"[A-Z]+\d+", r["linked_ids"]))
for r in csv.DictReader(open(BUILD / "accession_records.tsv"), delimiter="\t"):
    if r["accession"] in known:
        known.update(re.findall(r"(?:PRJ[A-Z]+|SRP|GSE)\d+", r["linked"] + " " + r["bioproject"] + " " + r["sra_link"]))
rows = []
for r in csv.DictReader(open(BUILD / "repo_sra_maize_cold.tsv"), delimiter="\t"):
    if r["bioproject"] not in known and r["srp"] not in known:
        rows.append({"accession": r["bioproject"], "kind": "SRA project", "date": r["create_date"][:10],
                     "n_runs": r["n_experiments"], "submitter": r["center"], "title": r["study_name"],
                     "example_samples": r["example_titles"]})
for r in csv.DictReader(open(BUILD / "repo_geo_maize_cold.tsv"), delimiter="\t"):
    if "high throughput" in r["gdstype"].lower() and r["accession"] not in known and r["bioproject"] not in known:
        rows.append({"accession": r["accession"], "kind": "GEO series", "date": r["pdat"][:10], "n_runs": r["n_samples"],
                     "submitter": "", "title": r["title"], "example_samples": f"PMID {r['pubmed']}" if r["pubmed"] else ""})
COLDISH = re.compile(r"cold|chill|low[- ]?temp|freez|frost|°C|temperature|highland", re.I)
NOT_RNA = re.compile(r"small RNA|miRNA|nucleosome|methylation|ChIP|ATAC|MOA-seq|CUT&Tag", re.I)
for r in rows:
    r["likely_cold_rnaseq"] = "yes" if COLDISH.search(r["title"] + " " + r["example_samples"]) and not NOT_RNA.search(r["title"]) else ""
rows.sort(key=lambda r: -int(r["n_runs"] or 0))
with open(BUILD / "uncatalogued_repository_records.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, ["accession", "kind", "date", "n_runs", "submitter", "title", "example_samples", "likely_cold_rnaseq"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"{len(rows)} unmatched records, {sum(1 for r in rows if r['likely_cold_rnaseq'])} likely cold RNA-seq")
