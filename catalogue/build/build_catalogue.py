#!/usr/bin/env python3
"""Merge the screening agents' extraction tables into the maize cold RNA-seq catalogue (M1).

Inputs : catalogue/build/agent_out/pubmed_batch_*.tsv, web_literature.tsv (agent extractions),
         manual_additions.tsv (studies added by hand from dataset leads),
         catalogue/build/agent_out/related_reports.tsv (RELATED items from the agents' reports)
Outputs: catalogue/maize_cold_rnaseq_studies.tsv  one row per included study
         catalogue/maize_cold_related.tsv         microarray / small RNA / reanalysis / ambiguous items
Usage (project root): python3 catalogue/build/build_catalogue.py
"""
import csv
import glob
import re
import sys
from collections import Counter
from pathlib import Path

BUILD = Path(__file__).resolve().parent
CAT = BUILD.parent

rows = []
for f in sorted(glob.glob(str(BUILD / "agent_out" / "pubmed_batch_*.tsv"))) + [str(BUILD / "agent_out" / "web_literature.tsv"), str(BUILD / "agent_out" / "manual_additions.tsv")]:
    if Path(f).exists():
        with open(f, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                r["_source_file"] = Path(f).name
                rows.append(r)


def key(d):
    pm = (d.get("pmid") or "").strip()
    return pm if pm.isdigit() else (d.get("doi_url") or d.get("title", "")).strip().lower()


# ---- manual curation (decisions made while reviewing the agents' reports) ----
# key (PMID or DOI URL) -> overrides: study_id, accession, acc_status, acc_note, fit ('core' /
# 'borderline' / 'tentative'), fit_note, tissue_group (list), system_group, related_to
OV = {
    # study_id clashes between batches: a/b in order of publication
    '39859301': dict(study_id='Sun2025a'), '40941777': dict(study_id='Sun2025b'),
    '36768951': dict(study_id='He2023a'),
    '34354054': dict(acc_note='Text cites PRJNA344653 (Zhang2017, reused?) and SRP327817 (probably its own deposit).'),
    '30342485': dict(acc_note='Paper text gives SRP162337; PRJNA492469 came from a web search.'), '37093290': dict(study_id='He2023b', fit='tentative',
                     fit_note='Abstract only; replicates, libraries and platform not confirmed.'),
    '35310673': dict(study_id='Meng2022a'),
    '35628306': dict(study_id='Meng2022b', fit='borderline',
                     fit_note='All 24 sequenced libraries are at 13 °C; the 25 °C comparator is not sequenced, so only time contrasts within cold.'),
    # datasets shared between papers (count the libraries once)
    '38794432': dict(related_to='Same BioProject (PRJNA1050059) as Dou2025.'),
    '41095130': dict(related_to='Same BioProject (PRJNA1050059) as Zhang2024; probably the same libraries.'),
    '39000567': dict(related_to='Same RNA-seq (N192/Ji853, 24 libraries) as Zhao2024b.'),
    '39273343': dict(related_to='Same RNA-seq as Zhao2024a; accession linked by design match.'),
    '33516847': dict(related_to='B144 vs Q319 at 5 °C, same group as Yu2023 and Yu2024.', acc_status='lead',
                     acc_note='PRJNA666026 came from a web search; the harvester finds it in Yu2023, not in this paper.'),
    '37628649': dict(related_to='B144 vs Q319, same group as Yu2021 and Yu2024; may reuse PRJNA666026.'),
    '39596336': dict(related_to='B144 vs Q319 at 5 °C 24 h, same group as Yu2021 and Yu2023; check for shared libraries.'),
    # fit
    '39766968': dict(fit='borderline', fit_note='Post-harvest storage of detached sweet corn kernels at 4 vs 28 °C, not cold stress of plants.'),
    '36616288': dict(fit='borderline', fit_note='Cold test per Methods: seeds soaked at 4 °C for 3 d, then germinated at 25 °C; check the design.'),
    '39095994': dict(fit='tentative', fit_note='Abstract only; a cold vs control RNA-seq design is not stated.'),
    '39701501': dict(fit='tentative', fit_note='Abstract only; temperatures, tissues and design not confirmed.'),
    '32505784': dict(fit='tentative', fit_note='Abstract only (AMF x low temperature); genotype, temperatures and design not confirmed.'),
    '31217070': dict(fit='tentative', fit_note='Abstract only; temperature, tissue and replicates not reported there.'),
    '41237517': dict(fit='tentative', fit_note='Abstract only; tissue, time points and libraries not confirmed.'),
    # found outside PubMed; only abstracts or snippets seen
    'https://doi.org/10.1016/j.cj.2017.03.005': dict(fit='tentative', fit_note='Not in PubMed; design from the abstract only.'),
    'https://doi.org/10.3390/agronomy15051067': dict(fit='tentative', fit_note='Not in PubMed; BSA-seq + transcriptome of RIL pools; design unverified.'),
    'https://doi.org/10.1101/2023.03.07.531599': dict(fit='tentative', fit_note='Preprint; design from search snippets only.'),
    'https://doi.org/10.64898/2026.03.04.709568': dict(fit='tentative', fit_note='Preprint; design from the abstract only.'),
}

TOK = re.compile(r"\b(PRJNA\d+|PRJEB\d+|PRJDB\d+|PRJCA\d+|SRP\d+|SRX\d+|SRR\d+|ERP\d+|DRP\d+|GSE\d+|CRA\d+|E-MTAB-\d+)\b")


def tissue_group(t):
    t = t.lower()
    g = []
    if re.search(r"root|radicle", t): g.append("Root")
    if re.search(r"leaf|leaves|blade|sheath", t): g.append("Leaf")
    if re.search(r"shoot|apex|meristem|mesocotyl|coleoptile|hypocotyl|stem", t): g.append("Shoot")
    if re.search(r"whole seedling|whole plant|^seedlings?\b(?! leaf| leaves| root)", t): g.append("Whole seedling")
    if re.search(r"\bseeds?\b|embryo|germinat|imbib|radicle", t): g.append("Seed/germination")
    if re.search(r"kernel|endosperm|ear|tassel|pollen|anther|silk", t) and "ear leaf" not in t: g.append("Reproductive")
    return g or ["NR"]


def system_group(s):
    s = s.lower()
    if "field" in s: return "Field"
    if re.search(r"greenhouse|glasshouse", s): return "Greenhouse"
    if re.search(r"chamber|incubator|climate|phytotron|controlled", s): return "Growth chamber"
    return "NR"


out, seen, dupes = [], {}, []
for d in rows:
    k = key(d)
    if k in seen:
        dupes.append((d.get("study_id", ""), seen[k]))
        continue
    seen[k] = d.get("study_id", "")
    o = dict(d)
    ov = OV.get(k, {})
    for f in ("study_id", "accession"):
        if f in ov:
            o[f] = ov[f]
    toks = list(dict.fromkeys(TOK.findall(o.get("accession", ""))))
    if "acc_status" in ov:
        st = ov["acc_status"]
    elif toks:
        st = "stated"
    elif "request" in o.get("accession", "").lower():
        st = "on request"
    else:
        st = "to resolve"
    src = o.get("source_of_details", "").lower()
    o.update(
        fit=ov.get("fit", "core"), fit_note=ov.get("fit_note", ""),
        details_level="abstract only" if src.startswith("abstract") and "pmc" not in src else "full text",
        tissue_group="; ".join(ov.get("tissue_group") or tissue_group(o.get("tissue", ""))),
        system_group=ov.get("system_group") or system_group(o.get("growth_system", "")),
        accession="; ".join(toks) or o.get("accession", ""),
        accession_status=st, accession_note=ov.get("acc_note", ""), related_to=ov.get("related_to", ""),
    )
    out.append(o)

# accessions harvested from the papers' PMC / Europe PMC full text and NCBI links
# (metadata/build/harvest_accessions.py -> accession_evidence.tsv); SRX/SRR run IDs left out
ev_path = CAT.parent / "metadata" / "build" / "accession_evidence.tsv"
harvest = {}
if ev_path.exists():
    for e in csv.DictReader(open(ev_path, encoding="utf-8"), delimiter="\t"):
        a_ = e["accession"]
        if re.match(r"(PRJ|GSE|SRP|ERP|DRP|CRA|E-MTAB)", a_):
            src = {"pmc_fulltext_scan": "PMC text", "epmc_fulltext_scan": "EPMC text"}.get(e["source"], e["source"].replace("elink_pubmed_", "NCBI link ").replace("epmc_annotation_", "EPMC annotation ").replace("epmc_dblink", "EPMC dblink"))
            harvest.setdefault(e["study_id"], {}).setdefault(a_, set()).add(src)
for o in out:
    h_ = harvest.get(o["study_id"], {})
    o["accession_harvested"] = "; ".join(f"{a_} ({', '.join(sorted(v))})" for a_, v in sorted(h_.items()))
    if o["accession_status"] == "to resolve" and h_:
        o["accession_status"] = "found by harvester"

ids = Counter(o["study_id"] for o in out)
clash = [i for i, n in ids.items() if n > 1]
if clash:
    sys.exit(f"duplicate study_id (set OV study_id): {clash}")
out.sort(key=lambda o: (-int(o["year"]) if str(o.get("year", "")).isdigit() else 0, o["study_id"]))

cols = ["study_id", "fit", "fit_note", "details_level", "year", "first_author", "journal", "title", "doi_url", "pmid", "pmcid",
        "design_category", "genotypes", "n_genotypes", "tissue", "tissue_group", "stage_age", "growth_system", "system_group",
        "cold_treatments", "control_temperature", "other_factors", "timepoints_duration", "bio_reps", "n_libraries",
        "platform_layout", "accession", "accession_status", "accession_harvested", "accession_note", "repository", "related_to", "other_omics",
        "key_finding", "notes", "source_of_details"]
with open(CAT / "maize_cold_rnaseq_studies.tsv", "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for o in out:
        w.writerow([str(o.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols])

rel = BUILD / "agent_out" / "related_reports.tsv"
if rel.exists():
    (CAT / "maize_cold_related.tsv").write_text(rel.read_text(encoding="utf-8"), encoding="utf-8")

print("studies:", len(out), "| duplicates merged:", len(dupes))
print("fit:", dict(Counter(o["fit"] for o in out)))
print("accession:", dict(Counter(o["accession_status"] for o in out)))
print("design:", dict(Counter(t.strip() for o in out for t in o.get("design_category", "").split("+") if t.strip())))
print("tissue:", dict(Counter(t for o in out for t in o["tissue_group"].split("; "))))
print("system:", dict(Counter(o["system_group"] for o in out)))
print("details:", dict(Counter(o["details_level"] for o in out)))
