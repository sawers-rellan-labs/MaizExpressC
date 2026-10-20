#!/usr/bin/env python3
"""Build metadata/runs.tsv: one row per sequencing run of every resolved dataset.

Inputs: metadata/accessions_resolved.tsv (resolved rows with role primary,
shared or partial), ENA run reports (metadata/build/ena_runs/), ENA sample
XML (submitter alias, title and attributes per BioSample), and for CNCB GSA
datasets the hand-captured metadata/build/gsa_runs.tsv plus gsa_md5_*.txt.

Include/exclude decisions are in RULES below, one function per study, so each
exclusion has a stated reason. Usage: python3 metadata/build/build_runs.py
"""
import csv
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harvest_accessions as h  # noqa: E402
import verify_accessions as v  # noqa: E402

ROOT = h.ROOT
BUILD = h.BUILD
OUT = ROOT / "metadata" / "runs.tsv"

COLS = ["run_accession", "study_id", "also_used_by", "bioproject", "biosample", "experiment",
        "sample_title", "organism", "library_strategy", "library_source", "library_layout",
        "read_length", "instrument", "n_spots", "n_bases", "fastq_urls", "fastq_md5",
        "biosample_attributes", "include", "exclude_reason"]
SKIP_TAGS = {"INSDC secondary accession", "INSDC center name", "INSDC first public",
             "INSDC last update", "INSDC status", "ENA-CHECKLIST", "ENA-FIRST-PUBLIC",
             "ENA-LAST-UPDATE", "SRA accession", "Submitter Id", "organism", "BioSampleModel"}


def sample_xml(accessions):
    """Return {sample accession: {alias, title, attributes}} from the ENA browser API."""
    out = {}
    accs = sorted(set(a for a in accessions if a))
    for i in range(0, len(accs), 50):
        txt = h.fetch("https://www.ebi.ac.uk/ena/browser/api/xml/" + ",".join(accs[i:i + 50]), "ebi")
        if not txt.strip():
            continue
        for s in ET.fromstring(txt).iter("SAMPLE"):
            attrs = {}
            for a in s.iter("SAMPLE_ATTRIBUTE"):
                tag, val = a.findtext("TAG", ""), a.findtext("VALUE", "")
                if tag not in SKIP_TAGS:
                    attrs[tag] = val
            out[s.get("accession")] = {"alias": s.get("alias", ""), "title": s.findtext("TITLE", ""),
                                       "attributes": attrs}
    return out


def ena_source(row):
    """ENA-queryable accession for a dataset row (GEO series go through their BioProject)."""
    acc = row["accession"]
    if acc.startswith(("PRJ", "SRP", "ERP", "DRP")):
        return acc
    m = re.search(r"PRJNA\d+", row["linked_ids"])
    return m.group(0) if m else acc


def load_ena(acc):
    path = BUILD / "ena_runs" / f"{acc}.tsv"
    if not path.exists():
        v.ena_runs(acc)
    return list(csv.DictReader(open(path), delimiter="\t")) if path.exists() else []


def read_length(r):
    try:
        n, b = int(r["read_count"]), int(r["base_count"])
    except ValueError:
        return ""
    if not n:
        return ""
    return str(round(b / n / (2 if r["library_layout"] == "PAIRED" else 1)))


# ---------------------------------------------------------------- rules
# Each rule gets the assembled row (with 'alias', 'title', 'attrs') and returns
# None to keep it, or an exclusion reason. Generic checks run first.
RNA_STRATEGIES = {"RNA-Seq"}


def generic(row):
    if row["organism"] != "Zea mays":
        return f"non-maize ({row['organism']})"
    if row["library_strategy"] not in RNA_STRATEGIES and row["study_id"] != "Ma2021a":
        return f"not RNA-seq ({row['library_strategy']})"
    return None


def mates_by_readcount(rows):
    """Map each read-2 run to its read-1 run, for datasets that uploaded R1 and R2 as
    separate runs with the same read count (Ma2021a: checked 17 Oct 2026, the earlier
    accession's reads are named n/1 and the later one's n/2, and the R2 reads are
    reverse-complement mates)."""
    seen = {}
    out = {}
    for r in sorted(rows, key=lambda r: r["run_accession"]):
        key = r["n_spots"]
        if key and key in seen:
            out[r["run_accession"]] = seen[key]
        else:
            seen[key] = r["run_accession"]
    return out


RULES = {}


def rule(study):
    def deco(f):
        RULES[study] = f
        return f
    return deco


def ena_library_names(acc):
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + acc
           + "&result=read_run&fields=run_accession,library_name&format=tsv")
    return {r["run_accession"]: r["library_name"]
            for r in csv.DictReader(h.fetch(url, "ebi").splitlines(), delimiter="\t")}


# Per-study run rules go here, one function per study, each returning an exclusion
# reason or None (examples in MaizExpressN: duplicate uploads, deferred data types,
# runs from treatments outside the project's scope):
#
# @rule("Study2024")
# def _(r):
#     if "heat" in r["sample_title"]:
#         return "heat treatment, no cold contrast"


# ---------------------------------------------------------------- build
def main(dump=False):
    datasets = [d for d in csv.DictReader(open(ROOT / "metadata" / "accessions_resolved.tsv"),
                                          delimiter="\t")
                if d["status"] == "resolved" and d["role"] in ("primary", "shared", "partial")]
    by_acc = {}
    for d in datasets:
        by_acc.setdefault(d["accession"], []).append(d)

    rows = []
    gsa = list(csv.DictReader(open(BUILD / "gsa_runs.tsv"), delimiter="\t"))
    md5 = {}
    for f in BUILD.glob("gsa_md5_*.txt"):
        for line in open(f):
            m, p = line.split()
            md5.setdefault(p.split("/")[2], []).append((p, m))

    for acc, ds in by_acc.items():
        owner = OWNER.get(acc, ds[0]["study_id"])
        others = ";".join(d["study_id"] for d in ds if d["study_id"] != owner)
        if ds[0]["repository"] == "CNCB GSA":
            for g in (x for x in gsa if x["cra"] == acc):
                files = sorted(md5.get(g["run_accession"], []))
                rows.append({
                    "run_accession": g["run_accession"], "study_id": owner, "also_used_by": others,
                    "bioproject": re.search(r"PRJCA\d+", ds[0]["linked_ids"]).group(0),
                    "biosample": g["sample_accession"], "experiment": "",
                    "sample_title": g["sample_title"], "organism": "Zea mays",
                    "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC",
                    "library_layout": g["layout"], "read_length": "", "instrument": g["instrument"],
                    "n_spots": "", "n_bases": "",
                    "fastq_urls": ";".join("https://download.cncb.ac.cn/" + g["ftp_dir"].split("/")[0]
                                           + p for p, _ in files),
                    "fastq_md5": ";".join(m for _, m in files),
                    "alias": g["run_alias"], "title": g["sample_title"],
                    "attrs": {"sample_name": g["sample_name"], "cultivar": g["cultivar"],
                              "tissue": g["tissue"], "repository": "CNCB GSA " + acc},
                })
            continue
        src = ena_source(ds[0])
        ena = load_ena(src)
        samples = sample_xml(r["sample_accession"] for r in ena)
        for r in ena:
            s = samples.get(r["sample_accession"], {"alias": "", "title": "", "attributes": {}})
            rows.append({
                "run_accession": r["run_accession"], "study_id": owner, "also_used_by": others,
                "bioproject": r["study_accession"], "biosample": r["sample_accession"],
                "experiment": r["experiment_accession"],
                "sample_title": r["sample_title"] or r["experiment_title"],
                "organism": r["scientific_name"], "library_strategy": r["library_strategy"],
                "library_source": r["library_source"], "library_layout": r["library_layout"],
                "read_length": read_length(r), "instrument": r["instrument_model"],
                "n_spots": r["read_count"], "n_bases": r["base_count"],
                "fastq_urls": ";".join("https://" + u for u in r["fastq_ftp"].split(";") if u),
                "fastq_md5": r["fastq_md5"],
                "alias": s["alias"], "title": s["title"], "attrs": s["attributes"],
            })

    for r in rows:
        split = SPLIT.get(r["bioproject"])
        if split:
            r["study_id"], r["also_used_by"] = split(r)
    apply_rules(rows)
    for r in rows:
        r["biosample_attributes"] = json.dumps(
            {"alias": r.pop("alias"), "title": r.pop("title"), **r.pop("attrs")}, ensure_ascii=False)
    tmp = OUT.with_suffix(".tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["study_id"], r["run_accession"])))
    tmp.replace(OUT)
    n_in = sum(r["include"] == "TRUE" for r in rows)
    print(f"{len(rows)} runs, {n_in} included -> {OUT}", file=sys.stderr)
    return rows


def apply_rules(rows):
    by_study = {}
    for r in rows:
        by_study.setdefault(r["study_id"], []).append(r)
    for study, rs in by_study.items():
        mates = MATES.get(study, lambda rs: {})(rs)
        by_acc = {r["run_accession"]: r for r in rs}
        for r2, r1 in mates.items():
            # the read-1 run carries both files, R1 first (make_manifests.py relies on the order)
            for col in ("fastq_urls", "fastq_md5"):
                by_acc[r1][col] = ";".join(x for x in (by_acc[r1][col], by_acc[r2][col]) if x)
        for r in rs:
            mate_of = mates.get(r["run_accession"])
            reason = (NOT_USED.get(study)
                      or generic(r)
                      or (f"read 2 of {mate_of} (R1 and R2 uploaded as separate runs); "
                          f"fetched as {mate_of}'s second file" if mate_of else None)
                      or (RULES.get(study, lambda r: None)(r)))
            r["include"] = "FALSE" if reason else "TRUE"
            r["exclude_reason"] = reason or ""
            extra = ALSO_USED.get(study)
            if extra and not reason:
                add = extra(r)
                if add:
                    r["also_used_by"] = ";".join(x for x in (r["also_used_by"], add) if x)


MATES = {}       # study -> function pairing R1/R2 uploaded as separate runs (see mates_by_readcount)

# Whole studies not used in the analysis: {study_id: "not used in this analysis"}.
# Reasons live in private/ (gitignored), never in committed files.
NOT_USED = {}

# Datasets claimed by two studies: who owns the runs.
OWNER = {}       # accession -> owning study_id
SPLIT = {}       # accession -> function(row) -> (owner, also_used_by)
ALSO_USED = {}   # owning study -> function(row) -> borrowing study or ""

if __name__ == "__main__":
    main()
