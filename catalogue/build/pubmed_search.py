#!/usr/bin/env python3
"""PubMed search for maize cold-stress transcriptome studies (M1).

Runs each query through NCBI E-utilities, records hit counts, and writes:
  catalogue/build/pubmed_queries.tsv     query number, query, hits, new PMIDs added
  catalogue/build/pubmed_candidates.tsv  union of PMIDs with title, journal, year, abstract
  catalogue/build/pubmed_prescreen_excluded.tsv  records set aside by the keyword pre-screen, with the reason
  catalogue/build/batches/batch_<n>.tsv  the remaining candidates split for screening

Responses are cached (metadata/build/cache), so reruns are offline.
Usage (project root): python3 catalogue/build/pubmed_search.py [batch_size]
"""
import csv
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "metadata" / "build"))
import harvest_accessions as h  # noqa: E402

OUT = Path(__file__).resolve().parent
MAIZE = '(maize OR "Zea mays" OR corn)'
COLD = '(cold OR chilling OR "low temperature" OR "low-temperature" OR freezing OR frost)'
TX = ('("RNA-seq" OR "RNA sequencing" OR "RNA-sequencing" OR transcriptome OR transcriptomic OR '
      'transcriptomics OR "differentially expressed" OR "gene expression profiling")')
QUERIES = [
    f"{MAIZE} AND {COLD} AND {TX}",
    '("Zea mays"[MeSH] OR maize) AND ("Cold Temperature"[MeSH] OR "Cold-Shock Response"[MeSH] OR '
    '"Freezing"[MeSH]) AND ("Gene Expression Profiling"[MeSH] OR "Transcriptome"[MeSH] OR '
    '"Sequence Analysis, RNA"[MeSH] OR "RNA-Seq"[MeSH] OR "Gene Expression Regulation, Plant"[MeSH])',
    f'{MAIZE} AND ("cold stress" OR "chilling stress" OR "cold tolerance" OR "chilling tolerance" OR '
    '"cold resistance" OR "chilling injury" OR "low temperature stress" OR "cold acclimation" OR '
    '"freezing tolerance" OR "cold-tolerant" OR "chilling-tolerant") AND (transcriptome OR transcriptomic '
    'OR "RNA-seq" OR "differentially expressed" OR "gene expression" OR "expression profiling")',
    f'{MAIZE} AND (cold OR chilling OR "low temperature") AND ("differentially expressed genes" OR DEGs '
    'OR WGCNA OR "co-expression" OR lncRNA OR "long non-coding" OR circRNA OR "single-cell" OR '
    '"single cell" OR "alternative splicing" OR eQTL OR TWAS)',
    f'{MAIZE} AND (cold OR chilling OR "low temperature") AND (germination OR emergence OR "seed vigor" '
    'OR "seed vigour" OR imbibition) AND (transcriptome OR transcriptomic OR "RNA-seq" OR '
    '"differentially expressed" OR "gene expression")',
    '(maize OR "Zea mays") AND ("abiotic stress" OR "abiotic stresses" OR "multiple stresses" OR '
    '"environmental stress" OR "environmental stresses") AND ("RNA-seq" OR transcriptome OR '
    'transcriptomic) AND (cold OR chilling OR "low temperature" OR temperature)',
    f'{MAIZE} AND ("suboptimal temperature" OR "sub-optimal temperature" OR "low night temperature" OR '
    '"night chilling" OR "early sowing" OR "early planting" OR "cool temperature" OR "cool temperatures" '
    'OR "temperature stress") AND (transcriptome OR transcriptomic OR "RNA-seq" OR '
    '"differentially expressed" OR "gene expression")',
    f'{MAIZE} AND (cold OR chilling OR "low temperature") AND (ICE1 OR CBF OR DREB OR DREB1 OR '
    '"cold-responsive" OR "cold responsive" OR "cold-regulated" OR "cold-induced") AND (transcriptome '
    'OR "RNA-seq" OR "differentially expressed" OR "gene expression profiling")',
    '(maize OR "Zea mays" OR teosinte OR Zea) AND (teosinte OR mexicana OR highland OR "high altitude" '
    'OR "high elevation" OR temperate OR tropical) AND (cold OR chilling OR "low temperature") AND '
    '(transcriptome OR transcriptomic OR "RNA-seq" OR "gene expression")',
    '(maize OR "Zea mays") AND (cold OR chilling OR "low temperature" OR freezing) AND (microRNA OR '
    'miRNA OR "small RNA" OR degradome OR microarray)',
]


def esearch(term):
    txt = h.eutil("esearch", db="pubmed", term=term, retmax=5000, retmode="json")
    res = json.loads(txt)["esearchresult"]
    return int(res["count"]), res.get("idlist", [])


def efetch(pmids):
    recs = {}
    for i in range(0, len(pmids), 150):
        xml = h.eutil("efetch", db="pubmed", id=",".join(pmids[i:i + 150]), retmode="xml")
        for art in ET.fromstring(xml).iter("PubmedArticle"):
            pmid = art.findtext(".//PMID")
            title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
            abst = " ".join("".join(a.itertext()) for a in art.iter("AbstractText"))
            year = (art.findtext(".//JournalIssue/PubDate/Year") or art.findtext(".//ArticleDate/Year")
                    or (art.findtext(".//JournalIssue/PubDate/MedlineDate") or "")[:4])
            doi = next((e.text for e in art.iter("ArticleId") if e.get("IdType") == "doi"), "")
            pmc = next((e.text for e in art.iter("ArticleId") if e.get("IdType") == "pmc"), "")
            ptypes = "; ".join(p.text for p in art.iter("PublicationType"))
            recs[pmid] = {"pmid": pmid, "year": year, "journal": art.findtext(".//Journal/Title") or "",
                          "title": re.sub(r"\s+", " ", title).strip(), "doi": doi, "pmcid": pmc,
                          "publication_types": ptypes, "abstract": re.sub(r"\s+", " ", abst).strip()}
    return recs


MAIZE_RE = re.compile(r"\b(maize|zea|corn|teosinte|mays)\b", re.I)
TEMP_RE = re.compile(r"(cold|chill|low[- ]temperature|freez|frost|sub-?optimal temperature|cool|°C|temperature)", re.I)
EXPR_RE = re.compile(r"(rna[- ]?seq|transcriptom|differentially expressed|\bDEGs?\b|sequencing|expression profil|"
                     r"microarray|mirna|small rna|degradome|gene expression|lncrna|co-expression|wgcna)", re.I)


def prescreen(r):
    if not r["abstract"]:
        return ""
    t = r["title"] + " " + r["abstract"]
    why = [m for m, rx in (("no Zea term", MAIZE_RE), ("no temperature term", TEMP_RE), ("no expression term", EXPR_RE))
           if not rx.search(t)]
    return "; ".join(why)


def main():
    size = int(sys.argv[1]) if len(sys.argv) > 1 else 36
    seen, rows = [], []
    for i, q in enumerate(QUERIES, 1):
        n, ids = esearch(q)
        new = [x for x in ids if x not in seen]
        seen += new
        rows.append({"query": i, "hits": n, "new": len(new), "text": q})
        print(f"query {i}: {n} hits, {len(new)} new")
    with open(OUT / "pubmed_queries.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, ["query", "hits", "new", "text"], delimiter="\t")
        w.writeheader(); w.writerows(rows)
    recs = efetch(seen)
    cols = ["pmid", "year", "journal", "title", "doi", "pmcid", "publication_types", "abstract"]
    allrecs = sorted(recs.values(), key=lambda r: (r["year"], r["pmid"]), reverse=True)
    with open(OUT / "pubmed_candidates.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t"); w.writeheader(); w.writerows(allrecs)
    # lenient keyword pre-screen: set aside records whose title + abstract lack a Zea term, any
    # temperature term or any expression term (records without an abstract go to screening)
    order, aside = [], []
    for r in allrecs:
        why = prescreen(r)
        (aside if why else order).append({**r, "reason": why} if why else r)
    with open(OUT / "pubmed_prescreen_excluded.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, ["pmid", "year", "title", "reason"], delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(aside)
    print(f"pre-screen: {len(order)} to screen, {len(aside)} set aside")
    (OUT / "batches").mkdir(exist_ok=True)
    for old in (OUT / "batches").glob("batch_*.tsv"):
        old.unlink()
    for b in range(0, len(order), size):
        with open(OUT / "batches" / f"batch_{b // size + 1}.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, cols, delimiter="\t"); w.writeheader(); w.writerows(order[b:b + size])
    print(f"{len(seen)} unique PMIDs, {len(recs)} records fetched, {-(-len(order) // size)} batches of up to {size}")


if __name__ == "__main__":
    main()
