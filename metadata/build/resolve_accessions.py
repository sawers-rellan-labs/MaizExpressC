#!/usr/bin/env python3
"""Write metadata/accessions_resolved.tsv (M1 checkpoint) and the unmatched repository records.

Inputs: metadata/build/accession_records.tsv (verify_accessions.py), the ENA run reports in
metadata/build/ena_runs/, and the decisions below (one entry per study and accession, made by
reading the papers' data statements and the repository records; evidence in `evidence`).
Status: resolved / ambiguous / unavailable / not_found.
Role: primary (the study's own data), shared (the study uses another study's data; libraries counted
once), partial (only some runs are maize cold RNA-seq), duplicate, exclude (not the study's data),
review (needs a decision from Rubén), not_used (study left out of the analysis by Rubén).
Usage (project root): python3 metadata/build/resolve_accessions.py
"""
import csv
from collections import Counter
from pathlib import Path

BUILD = Path(__file__).resolve().parent
META = BUILD.parent
ROOT = META.parent

P, S, PT, X, R = "primary", "shared", "partial", "exclude", "review"
D = [
    # study, accession, status, role, evidence_type, evidence, notes
    ("Avila2018", "SRP162337", "resolved", P, "paper_statement", "Paper text gives SRP162337 (Europe PMC / PMC full text)", "BioProject PRJNA492469. 24 runs."),
    ("Bozic2025", "PRJEB80094", "resolved", P, "paper_statement", "Data availability statement (PMC)", "ENA. 8 runs, one library per genotype x treatment x time (no replicates)."),
    ("Cao2022", "", "not_found", "", "paper_statement", "Data availability statement: 'Not applicable' (PMC web page)", "No deposit stated. Line-name search found large 18-599 projects (PRJNA957932, PRJNA644582) with no cold runs identified."),
    ("Dou2025", "PRJNA1050059", "resolved", S, "paper_statement", "Paper states PRJNA1050059 (PMC)", "Same BioProject as Zhang2024 (12 runs); libraries owned by Zhang2024."),
    ("Du2023", "", "not_found", "", "paper_statement", "Data availability statement: 'Not applicable' (PMC web page)", "No deposit stated; no repository match for Ye478 low-temperature germination."),
    ("Frey2020", "PRJNA556806", "resolved", P, "paper_statement", "Paper and Europe PMC annotation", "176 runs: 21 DH lines plus reference lines PE0075 and B104, cold vs control root tips."),
    ("Friero2023", "PRJNA877073", "ambiguous", R, "repository_match", "BioProject title 'Molecular signatures of temperature-dependent root meristem growth in maize' matches the paper's design; publisher page not readable (ScienceDirect)", "14 runs (DR_0-30, DR_12-20, DR_12-30 x 3 reps and more). Confirm from the paper's data statement."),
    ("Gu2017", "PRJNA361019", "ambiguous", "not_used", "harvested", "SRX2672484 found by the accession harvester; it belongs to PRJNA361019 (SRP102470)", "Only 1 run deposited (He344 cold leaf, HiSeq 4000) against the paper's 6 libraries. Not used in this analysis (Rubén, 20 Oct 2026)."),
    ("He2016", "SRP093477", "resolved", P, "paper_statement", "Paper (PMC)", "5 runs; no biological replication."),
    ("He2023a", "CRA007319", "resolved", P, "repository_pubmed_link", "CNCB GSA record links PMID 36768951 (checked in a browser)", "CNCB GSA, BioProject PRJCA010234. 36 runs (72 files, 107 GB): B73 and introgression line IB030, control and cold time points."),
    ("Han2026", "CRA027029", "unavailable", P, "paper_statement", "Data availability statement on the Springer page: GSA CRA027029", "Not released in GSA (checked 20 Oct 2026). PRJNA1314463 (a repository match tried earlier) is not this study."),
    ("Jia2020", "PRJNA598312", "resolved", PT, "repository_match", "SRA titles 'RNA-seq of Zea mays: wild type leaf ... at 28' etc. (University of Nebraska, RDM4); the paper's 18 leaf libraries are WT and rdm4 x 28 / 10 / 2 °C x 3 (accepted manuscript)", "26 runs: 18 leaf (cold series), 6 endosperm (no cold), 2 exome capture. Keep the 18 leaf runs."),
    ("Jia2020", "PRJNA306879", "resolved", R, "paper_statement", "Accession cited in the paper (publisher PDF via subagent); Europe PMC", "34 runs: 12 RNA-seq (probably the mapping data), 22 exome. The cold leaf series is in PRJNA598312; check at M2 whether any of these 12 belong to it."),
    ("Jiang2023", "PRJNA911416", "resolved", P, "harvested", "PMC full-text scan", "18 runs (BioProject titled 'Genome sequencing and assembly' but all runs are RNA-Seq)."),
    ("Li2016", "SRP065916", "resolved", P, "paper_statement", "Europe PMC annotation", "12 runs (KR701 and Hei8834, control vs freezing)."),
    ("Li2017", "SRP080208", "resolved", P, "paper_statement", "Paper (PMC web page)", "12 runs."),
    ("Li2020a", "PRJNA645274", "resolved", PT, "paper_statement", "Paper (PMC)", "42 runs: cold (4/10/16 °C vs 25 °C) and heat series; keep the cold series and controls."),
    ("Li2021", "PRJNA533857", "unavailable", P, "paper_statement", "Paper states PRJNA533857 (PMC)", "BioProject exists but has no released runs in SRA or ENA (checked 20 Oct 2026)."),
    ("Li2022a", "", "not_found", "", "paper_statement", "Data availability statement: 'No New data' (PMC web page)", "RNA-seq in the paper may reuse earlier data; source not identified."),
    ("Li2022b", "PRJNA822071", "resolved", "not_used", "harvested", "PMC full-text scan", "12 runs, Oxford Nanopore MinION full-length reads (single). Not used in this analysis (Rubén, 20 Oct 2026)."),
    ("Li2022c", "PRJNA778546", "resolved", P, "paper_statement", "Paper (PMC web page)", "15 runs."),
    ("Li2022d", "PRJNA759333", "resolved", P, "harvested", "PMC full-text scan", "24 runs (two near-isogenic lines)."),
    ("Li2025", "GSE299259", "resolved", P, "harvested", "PMC full-text scan", "18 runs (zmzfp69 mutant vs WT)."),
    ("Liu2026", "PRJNA1393347", "resolved", P, "harvested", "PMC full-text scan", "84 runs, matching the paper's 84 libraries (cold, heat, drought)."),
    ("Lu2017", "GSE76939", "resolved", P, "paper_statement", "Paper (PMC) and NCBI PubMed link", "Teosinte (Zea mays subsp. mexicana), 6 runs; BioProject PRJNA309111, SRP068562."),
    ("Lu2017", "PRJNA866469", "resolved", X, "agent", "Listed by the screening agent", "Zea luxurians genome project; not this study's RNA-seq."),
    ("Makarevitch2015", "PRJNA244661", "resolved", P, "paper_statement", "Data availability statement (PMC)", "33 runs (cold, heat, salt, UV; B73, Mo17, Oh43 and later B37/PH207 runs used by Waters2017). Corrected S1 table (erratum PMC4599801) maps runs to samples."),
    ("Meng2022a", "PRJNA795360", "resolved", P, "paper_statement", "Paper (PMC)", "12 runs."),
    ("Meng2025", "PRJNA1088921", "resolved", P, "paper_statement", "Data availability statement (Wiley page)", "12 runs."),
    ("Meng2025", "PRJNA1088971", "resolved", P, "paper_statement", "Data availability statement (Wiley page)", "12 runs."),
    ("Meng2025", "PRJNA1126396", "resolved", P, "paper_statement", "Data availability statement (Wiley page)", "12 runs."),
    ("Meng2025", "PRJNA1126430", "resolved", P, "paper_statement", "Data availability statement (Wiley page)", "12 runs."),
    ("Qin2025", "", "not_found", "", "paper_statement", "Data availability statement: data included in the article (PMC web page)", "No deposit stated."),
    ("Shi2025", "", "not_found", "", "paper_statement", "Data availability: measurements in the supplemental file (PMC web page)", "No deposit stated."),
    ("Sommer2025", "PRJNA1243569", "resolved", P, "paper_statement", "Paper (PMC), Europe PMC annotation", "15 runs (B73 root hairs, control / mild / severe cold x 5 reps)."),
    ("Sun2025a", "", "not_found", "", "paper_statement", "Data availability statement: data contained within the article (PMC web page)", "No deposit stated."),
    ("Sun2025b", "", "not_found", "", "paper_statement", "Data availability statement: data in the article and supplement (PMC web page)", "No deposit stated."),
    ("Tang2026", "CRA041093", "unavailable", P, "harvested", "PMC full-text scan", "Not released in CNCB GSA (checked 20 Oct 2026)."),
    ("Tang2026", "CRA041107", "unavailable", P, "harvested", "PMC full-text scan", "Not released in CNCB GSA (checked 20 Oct 2026)."),
    ("Wang2023", "PRJNA767196", "resolved", PT, "paper_statement", "Data descriptor (Scientific Data); NCBI SRA", "132 runs from 11 species; 12 Zea mays runs (4 °C at 0/2/24/168 h x 3)."),
    ("Waters2017", "PRJNA244661", "resolved", S, "paper_statement", "Accession numbers section (Wiley page): PRJNA244661", "Same BioProject as Makarevitch2015; Waters2017 adds genotypes (B37, PH207). Runs are curated once under Makarevitch2015 with also_used_by."),
    ("Waters2017", "PRJNA385873", "ambiguous", R, "repository_match", "University of Minnesota project with inbred x hybrid cold runs; not cited in the paper", "54 runs (e.g. Oh43xB73, B73xPH207 cold reps). May belong to a later paper; listed for triage."),
    ("Xue2021", "GSE167156", "resolved", P, "repository_pubmed_link", "GEO record links PMID 34240135", "16 runs: four lines, chilling (BioProject PRJNA703398)."),
    ("Xue2021", "GSE167881", "resolved", P, "repository_pubmed_link", "GEO record links PMID 34240135", "48 runs: tissue-level local vs distal chilling (BioProject PRJNA705456)."),
    ("Yang2025", "PRJNA1111354", "resolved", P, "harvested", "PMC full-text scan", "15 runs (ABA/BR x cold)."),
    ("Yu2021", "PRJNA666026", "resolved", P, "repository_match", "B144/Q319 BioProject registered Sep 2020; 18 runs fit Yu2021's design (2 lines x 0/12/24 h x 3); Yu2023 cites it", "Owner of the B144/Q319 libraries; Yu2023 and Yu2024 reuse them."),
    ("Yu2023", "PRJNA666026", "resolved", S, "paper_statement", "PMC full-text scan and Europe PMC annotation", "Reuses Yu2021's libraries."),
    ("Yu2024", "PRJNA666026", "resolved", S, "paper_statement", "Data availability: 'No new data were created' (PMC web page); B144/Q319 project", "Reuses Yu2021's libraries."),
    ("Zeng2021", "SRP327817", "resolved", P, "harvested", "PMC full-text scan", "8 runs (ZmRR1-dependent cold responses)."),
    ("Zeng2021", "PRJNA344653", "resolved", X, "harvested", "Cited in the text", "Zhang2017's dataset, cited for comparison; not Zeng2021's data."),
    ("Zhang2017", "PRJNA344653", "resolved", PT, "repository_pubmed_link", "NCBI PubMed link", "108 runs: maize, sorghum and setaria; 36 Zea mays runs (cold time course)."),
    ("Zhang2020", "GSE146666", "resolved", P, "paper_statement", "GEO record links the PMID; Europe PMC annotation", "12 runs (BioProject PRJNA611589)."),
    ("Zhang2024", "PRJNA1050059", "resolved", P, "paper_statement", "Paper (PMC)", "12 runs; also cited by Dou2025."),
    ("Zhao2021", "PRJNA686250", "resolved", P, "paper_statement", "Data availability statement lists SRR13312968-72 (PMC web page); the runs belong to PRJNA686250", "6 runs (He344 roots, cold vs control)."),
    ("Zhao2024a", "PRJNA1096968", "unavailable", P, "paper_statement", "Paper (PMC)", "Not public in NCBI (checked 20 Oct 2026)."),
    ("Zhao2024a", "PRJNA1121054", "unavailable", P, "paper_statement", "Paper (PMC)", "Not public in NCBI (checked 20 Oct 2026)."),
    ("Zhao2024a", "PRJNA917698", "resolved", X, "paper_statement", "Cited in the paper", "Light-quality mesocotyl study (24 runs), not cold."),
    ("Zhao2024b", "PRJNA1096968", "unavailable", S, "agent", "Linked by design match to Zhao2024a", "Same RNA-seq as Zhao2024a; not public."),
    ("Zhou2022", "PRJNA747925", "resolved", PT, "paper_statement", "Paper (PMC)", "291 runs: heat and cold time courses across genotypes; keep cold and controls."),
    ("Zhou2024", "PRJNA852372", "resolved", P, "paper_statement", "Accession numbers section (OUP page)", "160 runs, DNBSEQ: cell-type / tissue-resolved root samples (cortex, endodermis, epidermis, pericycle) under cold."),
    ("Zhou2026", "PRJNA1378803", "unavailable", P, "paper_statement", "Paper (PMC, Europe PMC)", "Not public in NCBI (checked 20 Oct 2026)."),
]


def maize_rnaseq_runs(acc):
    """Maize RNA-seq runs in the ENA run report (organism unset counts as maize for maize-only projects)."""
    f = BUILD / "ena_runs" / f"{acc}.tsv"
    if not acc or not f.exists():
        return ""
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    return str(sum(1 for r in rows if r.get("library_strategy") == "RNA-Seq"
                   and (r.get("scientific_name", "").startswith("Zea") or not r.get("scientific_name"))))


def main():
    recs = {r["accession"]: r for r in csv.DictReader(open(BUILD / "accession_records.tsv"), delimiter="\t")}
    cat = [r["study_id"] for r in csv.DictReader(open(ROOT / "catalogue" / "maize_cold_rnaseq_studies.tsv"), delimiter="\t")]
    missing = sorted(set(cat) - {d[0] for d in D})
    if missing:
        raise SystemExit(f"no decision for: {missing}")
    cols = ["study_id", "accession", "accession_type", "repository", "linked_ids", "status", "role", "evidence_type",
            "n_maize_rnaseq_runs", "evidence", "notes"]
    out = []
    for sid, acc, st, role, et, ev, notes in D:
        r = recs.get(acc, {})
        rep = ("CNCB GSA" if acc.startswith(("CRA", "PRJCA")) else "ENA" if acc.startswith(("PRJEB", "ERP"))
               else "NCBI GEO" if acc.startswith("GSE") else "NCBI SRA" if acc else "")
        n = maize_rnaseq_runs(acc)
        if not n and r.get("n_runs") and r.get("organisms", "").startswith("Zea"):
            n = r["n_runs"]   # GEO series: ENA report filed under the linked SRA study
        if acc == "CRA007319":
            n = "36"   # CNCB GSA page (72 files, paired), checked in a browser
        out.append({"study_id": sid, "accession": acc, "accession_type": r.get("type", ""), "repository": rep,
                    "linked_ids": r.get("linked", ""), "status": st, "role": role, "evidence_type": et,
                    "n_maize_rnaseq_runs": n, "evidence": ev, "notes": notes})
    with open(META / "accessions_resolved.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t")
        w.writeheader()
        w.writerows(out)
    per_study = {}
    order = ["resolved", "ambiguous", "unavailable", "not_found"]
    for o in out:
        if o["role"] in ("exclude", "not_used"):
            continue
        s = per_study.get(o["study_id"])
        per_study[o["study_id"]] = o["status"] if s is None else min(s, o["status"], key=order.index)
    print("rows:", len(out), "| per study:", dict(Counter(per_study.values())))
    print("roles:", dict(Counter(o["role"] for o in out if o["role"])))


if __name__ == "__main__":
    main()
