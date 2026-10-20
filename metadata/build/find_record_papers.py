#!/usr/bin/env python3
"""Find papers for the unmatched maize cold repository records (likely_cold_rnaseq = yes).

For each record: its linked IDs (BioProject, SRA study, GEO series) from the ENA run report / NCBI,
then a Europe PMC full-text search for each ID, plus the NCBI BioProject -> PubMed link.
Writes metadata/build/record_papers.tsv (accession, ids searched, hits: PMID/PMCID/DOI, title, year).
Usage (project root): python3 metadata/build/find_record_papers.py
"""
import csv
import io
import json
import sys
import urllib.parse
from pathlib import Path

BUILD = Path(__file__).resolve().parent
sys.path.insert(0, str(BUILD))
import harvest_accessions as h  # noqa: E402

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"


def linked_ids(acc):
    ids = {acc}
    if acc.startswith(("PRJ", "SRP", "ERP", "DRP")):
        t = h.fetch(f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run"
                    f"&fields=study_accession,secondary_study_accession&format=tsv&limit=5", "ebi")
        for r in csv.DictReader(io.StringIO(t or ""), delimiter="\t"):
            ids.update(x for x in (r.get("study_accession"), r.get("secondary_study_accession")) if x)
    return sorted(ids)


def epmc(q):
    t = h.fetch(EPMC + urllib.parse.urlencode({"query": f'"{q}"', "format": "json", "resultType": "lite", "pageSize": 25}), "ebi")
    if not t:
        return []
    return [(r.get("pmid", ""), r.get("pmcid", ""), r.get("doi", ""), r.get("title", ""), r.get("pubYear", ""), r.get("source", ""))
            for r in json.loads(t).get("resultList", {}).get("result", [])]


def bioproject_pubmed(acc):
    if not acc.startswith("PRJNA"):
        return []
    txt = h.eutil("esearch", db="bioproject", term=f"{acc}[Project Accession]", retmode="json")
    ids = json.loads(txt)["esearchresult"].get("idlist", []) if txt else []
    out = []
    for u in ids:
        x = h.eutil("elink", dbfrom="bioproject", db="pubmed", id=u, retmode="json")
        for ls in (json.loads(x).get("linksets", []) if x else []):
            for d in ls.get("linksetdbs", []):
                out += [str(p) for p in d.get("links", [])]
    return out


rows = [r for r in csv.DictReader(open(BUILD / "uncatalogued_repository_records.tsv"), delimiter="\t") if r["likely_cold_rnaseq"] == "yes"]
out = []
for r in rows:
    ids = linked_ids(r["accession"])
    hits = {}
    for i in ids:
        for pm, pmc, doi, title, year, src in epmc(i):
            hits[pm or doi or title] = (pm, pmc, doi, title, year, f"Europe PMC full text cites {i}")
    for pm in bioproject_pubmed(r["accession"]):
        hits.setdefault(pm, (pm, "", "", "", "", "NCBI BioProject -> PubMed link"))
    if not hits:
        out.append({"accession": r["accession"], "record_title": r["title"], "ids_searched": "; ".join(ids), "pmid": "", "pmcid": "",
                    "doi": "", "paper_title": "", "year": "", "how_found": "none"})
    for pm, pmc, doi, title, year, how in hits.values():
        out.append({"accession": r["accession"], "record_title": r["title"], "ids_searched": "; ".join(ids), "pmid": pm, "pmcid": pmc,
                    "doi": doi, "paper_title": title, "year": year, "how_found": how})
with open(BUILD / "record_papers.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, list(out[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(out)
print(f"{len(rows)} records -> {sum(1 for o in out if o['pmid'] or o['doi'])} paper hits; "
      f"{len({o['accession'] for o in out if o['how_found'] == 'none'})} records with none")
