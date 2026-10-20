#!/usr/bin/env python3
"""Collect candidate data accessions per catalogue study (M1, accession resolution).

Sources, per study:
  - the catalogue's `accession` column (what the screening agents read in the paper or found by web search);
  - metadata/build/accession_evidence.tsv (harvest_accessions.py: PMC / Europe PMC full-text scans,
    Europe PMC annotations, NCBI PubMed links);
  - EXTRA below (accessions found by hand, e.g. title matches in the repository search tables).
SRX/SRR run IDs are left out when a study-level ID exists.
Writes metadata/build/accession_candidates.tsv (study_id, accession, source) for verify_accessions.py.
Usage (project root): python3 metadata/build/build_candidates.py
"""
import csv
import re
from collections import defaultdict
from pathlib import Path

BUILD = Path(__file__).resolve().parent
ROOT = BUILD.parents[1]
TOK = re.compile(r"\b(PRJNA\d+|PRJEB\d+|PRJDB\d+|PRJCA\d+|SRP\d+|ERP\d+|DRP\d+|SRX\d+|GSE\d+|CRA\d+|E-MTAB-\d+|E-MEXP-\d+)\b")
SRC = {"pmc_fulltext_scan": "PMC full text", "epmc_fulltext_scan": "Europe PMC full text",
       "epmc_annotation_MED": "Europe PMC annotation", "epmc_annotation_PMC": "Europe PMC annotation",
       "epmc_dblink": "Europe PMC database link", "elink_pubmed_bioproject": "NCBI PubMed link",
       "elink_pubmed_sra": "NCBI PubMed link", "elink_pubmed_gds": "NCBI PubMed link"}

# study_id -> [(accession, source)] found by hand
EXTRA = {
    "Zhao2021": [("PRJNA686250", "paper statement lists SRR13312968-72 (PMC web page); runs belong to PRJNA686250")],
    "Waters2017": [("PRJNA385873", "repository match: University of Minnesota, inbred x hybrid cold/heat runs (54)"),
                   ("PRJNA244661", "repository match: same lab's B37/Oh43/B73/Mo17 cold runs")],
    "Friero2023": [("PRJNA877073", "repository match: 'Molecular signatures of temperature-dependent root meristem growth in maize'")],
    "Meng2025": [("PRJNA1088971", "repository match: Jilin Agricultural University, maize seed cold + high humidity"),
                 ("PRJNA1088921", "repository match: Jilin Agricultural University, maize seed cold + high humidity"),
                 ("PRJNA1126396", "repository match: Jilin Agricultural University, maize seed cold + high humidity"),
                 ("PRJNA1126430", "repository match: Jilin Agricultural University, maize seed cold + high humidity")],
    "Han2026": [("CRA027029", "data availability statement on the Springer page")],
    "Zhou2024": [("PRJNA852372", "accession numbers section on the OUP page")],
    "Gu2017": [("PRJNA361019", "BioProject of SRX2672484 (found by the harvester)")],
    "Jia2020": [("PRJNA598312", "University of Nebraska RDM4 project; leaf RNA-seq at 28, 10 and 2 °C (SRA titles), matching the paper's design")],
    "Yu2024": [("PRJNA666026", "repository match: B144/Q319 project; the paper says no new data were created")],
}


def main():
    cat = list(csv.DictReader(open(ROOT / "catalogue" / "maize_cold_rnaseq_studies.tsv"), delimiter="\t"))
    ids = {r["study_id"] for r in cat}
    cand = defaultdict(lambda: defaultdict(set))
    for r in cat:
        src = "catalogue (agent: paper)" if r["accession_status"] == "stated" else "catalogue (agent: web search, unverified)"
        if "web" in r["notes"].lower() and r["accession_status"] == "stated" and "web" in r["source_of_details"].lower():
            src = "catalogue (agent: paper or web search)"
        for a in TOK.findall(r["accession"]):
            cand[r["study_id"]][a].add(src)
    for e in csv.DictReader(open(BUILD / "accession_evidence.tsv"), delimiter="\t"):
        if e["study_id"] in ids and TOK.fullmatch(e["accession"]):
            cand[e["study_id"]][e["accession"]].add(SRC.get(e["source"], e["source"]))
    for sid, pairs in EXTRA.items():
        for a, s in pairs:
            cand[sid][a].add(s)
    rows = []
    for sid in sorted(cand):
        accs = cand[sid]
        study_level = [a for a in accs if not a.startswith("SRX")]
        for a in sorted(accs):
            if a.startswith("SRX") and study_level:
                continue
            rows.append({"study_id": sid, "accession": a, "source": "; ".join(sorted(accs[a]))})
    with open(BUILD / "accession_candidates.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, ["study_id", "accession", "source"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    none = sorted(ids - set(cand))
    print(f"{len(rows)} candidates for {len(cand)} studies; no candidate: {len(none)}: {', '.join(none)}")


if __name__ == "__main__":
    main()
