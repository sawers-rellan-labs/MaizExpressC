#!/usr/bin/env python3
"""Find papers through datasets (M1): link maize cold BioProjects, GEO series and dataset leads to PubMed.

Some key studies never say "cold" or "RNA-seq" in their abstracts (e.g. multi-stress papers that
only say "abiotic stress"), so keyword searches miss them. This script takes every BioProject and
GEO series from metadata/build/repo_*_maize_cold.tsv plus the lead accessions in
agent_out/to_check.tsv, links them to PubMed with NCBI elink, and writes the papers not already
among the PubMed candidates to batches/batch_dataset.tsv for screening.
Usage (project root): python3 catalogue/build/dataset_papers.py
"""
import csv
import json
import re
import sys
from pathlib import Path

BUILD = Path(__file__).resolve().parent
ROOT = BUILD.parents[1]
sys.path.insert(0, str(ROOT / "metadata" / "build"))
import harvest_accessions as h  # noqa: E402
sys.path.insert(0, str(BUILD))
import pubmed_search as ps  # noqa: E402

MB = ROOT / "metadata" / "build"
projects, geo_pmids, leads = set(), {}, set()
for r in csv.DictReader(open(MB / "repo_bioproject_maize_cold.tsv"), delimiter="\t"):
    projects.add(r["accession"])
for r in csv.DictReader(open(MB / "repo_sra_maize_cold.tsv"), delimiter="\t"):
    if r["bioproject"].startswith("PRJ"):
        projects.add(r["bioproject"])
for r in csv.DictReader(open(MB / "repo_geo_maize_cold.tsv"), delimiter="\t"):
    for p in filter(None, r["pubmed"].split(";")):
        geo_pmids.setdefault(p, []).append(r["accession"])
    if r["bioproject"]:
        projects.add(r["bioproject"])
tc = BUILD / "agent_out" / "to_check.tsv"
if tc.exists():
    for r in csv.DictReader(open(tc), delimiter="\t"):
        leads.update(re.findall(r"\b(PRJ[NED][A-Z]?\d+|SRP\d+|SRX\d+|GSE\d+)\b", r["label"]))


def uids(db, term):
    txt = h.eutil("esearch", db=db, term=term, retmode="json")
    return json.loads(txt)["esearchresult"].get("idlist", []) if txt else []


def to_pubmed(acc):
    """PMIDs linked to a BioProject, SRA or GEO accession."""
    if acc.startswith("PRJ"):
        db, ids = "bioproject", uids("bioproject", f"{acc}[Project Accession]")
    elif acc.startswith(("SRP", "SRX")):
        db, ids = "sra", uids("sra", acc)
    else:
        db, ids = "gds", uids("gds", f"{acc}[ACCN] AND gse[ETYP]")
    out = set()
    for u in ids[:20]:
        txt = h.eutil("elink", dbfrom=db, db="pubmed", id=u, retmode="json")
        for ls in (json.loads(txt).get("linksets", []) if txt else []):
            for d in ls.get("linksetdbs", []):
                out.update(map(str, d.get("links", [])))
    return out


found = {}
for acc in sorted(projects | leads):
    for p in to_pubmed(acc):
        found.setdefault(p, set()).add(acc)
for p, accs in geo_pmids.items():
    found.setdefault(p, set()).update(accs)

screened = {r["pmid"] for r in csv.DictReader(open(BUILD / "pubmed_candidates.tsv"), delimiter="\t")}
new = sorted(p for p in found if p not in screened)
recs = ps.efetch(new) if new else {}
cols = ["pmid", "year", "journal", "title", "doi", "pmcid", "publication_types", "abstract", "linked_datasets"]
with open(BUILD / "batches" / "batch_dataset.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, cols, delimiter="\t")
    w.writeheader()
    for p in new:
        if p in recs:
            w.writerow({**recs[p], "linked_datasets": "; ".join(sorted(found[p]))})
with open(BUILD / "dataset_pubmed_links.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["pmid", "datasets", "already_screened"])
    for p in sorted(found):
        w.writerow([p, "; ".join(sorted(found[p])), p in screened])
print(f"{len(projects)} BioProjects + {len(leads)} leads + {len(geo_pmids)} GEO-linked papers -> "
      f"{len(found)} papers, {len(new)} not yet screened -> batches/batch_dataset.tsv")
