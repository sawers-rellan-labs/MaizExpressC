#!/usr/bin/env python3
"""Write per-study download manifests and nf-core/rnaseq samplesheets.

For each study given (or --slice for the M3 vertical slice), from
metadata/runs.tsv (include = TRUE) and metadata/samples.tsv (include = TRUE):

  pipeline/manifests/<study>.downloads.tsv  run, file_url, md5, dest (one row per FASTQ file)
  pipeline/manifests/<study>.samplesheet.csv  nf-core/rnaseq: sample,fastq_1,fastq_2,strandedness

Destination paths are on Hazel scratch (FASTQ_ROOT). Several runs of one
sample get one samplesheet row each with the same sample name; nf-core/rnaseq
concatenates them. Strandedness is 'auto' (inferred by salmon).

Usage: python3 pipeline/make_manifests.py <study> [<study> ...]   (or --all)
Single-file paired runs of SPLIT_HALVES studies also get pipeline/manifests/<study>.split.tsv.
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "pipeline" / "manifests"
FASTQ_ROOT = "/share/maize/rrellan/MaizeExpressC/fastq"
SLICE = []
# paired runs that ENA serves as one file holding all R1 then all R2 (see split_halves.sh)
SPLIT_HALVES = {"Zhou2024", "Friero2023"}


def pick_files(run):
    """Return [(url, md5), ...] for one run: one file (single) or R1, R2 (paired)."""
    urls = [u for u in run["fastq_urls"].split(";") if u]
    md5s = run["fastq_md5"].split(";")
    files = list(zip(urls, md5s + [""] * (len(urls) - len(md5s))))
    r1 = [f for f in files if re.search(r"_(1|f1|r1)\.f(ast)?q\.gz$", f[0])]
    r2 = [f for f in files if re.search(r"_(2|r2)\.f(ast)?q\.gz$", f[0])]
    if len(r1) == 1 and len(r2) == 1:
        return r1 + r2   # paired files, whatever the declared layout (Zhao2021 says SINGLE)
    if run["library_layout"] == "PAIRED":
        # ENA: _1/_2.fastq.gz; CNCB GSA: _f1/_r2 or _r1/_r2, .fq.gz or .fastq.gz
        r1 = [f for f in files if re.search(r"_(1|f1|r1)\.f(ast)?q\.gz$", f[0])]
        r2 = [f for f in files if re.search(r"_(2|r2)\.f(ast)?q\.gz$", f[0])]
        if len(r1) == 1 and len(r2) == 1:
            return r1 + r2
        if len(files) == 2 and not r1 and not r2:
            # Ma2021a: R1 and R2 were uploaded as separate runs; build_runs.py lists
            # the mate run's file second
            return files
        if len(files) == 1 and run["study_id"] in SPLIT_HALVES:
            return files   # one file: split into R1 / R2 after download
        if len(files) == 1:
            raise ValueError(f"{run['run_accession']}: paired run served by ENA as one file")
        raise ValueError(f"{run['run_accession']}: cannot find R1/R2 among {urls}")
    if len(files) == 1:
        return files
    raise ValueError(f"{run['run_accession']}: single-end run with {len(files)} files")


def main(studies):
    runs = {r["run_accession"]: r for r in csv.DictReader(
        open(ROOT / "metadata" / "runs.tsv"), delimiter="\t") if r["include"] == "TRUE"}
    run_sample = {}
    for s in csv.DictReader(open(ROOT / "metadata" / "samples.tsv"), delimiter="\t"):
        if s["include"] == "TRUE":
            for r in s["run_accessions"].split(";"):
                run_sample[r] = s["sample_id"]
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for study in studies:
        rs = sorted((r for r in runs.values() if r["study_id"] == study and
                     r["run_accession"] in run_sample), key=lambda r: r["run_accession"])
        if not rs:
            sys.exit(f"{study}: no included runs")
        dl, sheet, layouts, split = [], [], defaultdict(set), []
        total = 0
        for r in rs:
            files = pick_files(r)
            dest = [f"{FASTQ_ROOT}/{study}/{u.rsplit('/', 1)[1]}" for u, _ in files]
            for (u, m), d in zip(files, dest):
                dl.append([r["run_accession"], u, m, d])
            sample = run_sample[r["run_accession"]]
            if len(files) == 1 and r["library_layout"] == "PAIRED" and study in SPLIT_HALVES:
                base = dest[0].rsplit(".fastq.gz", 1)[0]
                split.append([r["run_accession"], dest[0], r["n_spots"], base + "_1.fastq.gz", base + "_2.fastq.gz"])
                dest = [base + "_1.fastq.gz", base + "_2.fastq.gz"]
            layouts[sample].add("PAIRED" if len(dest) == 2 else "SINGLE")
            sheet.append([sample, dest[0], dest[1] if len(dest) > 1 else "", "auto"])
            total += int(r["n_bases"] or 0)
        mixed = [s for s, l in layouts.items() if len(l) > 1]
        if mixed:
            sys.exit(f"{study}: samples mix single and paired runs: {mixed}")
        with open(OUTDIR / f"{study}.downloads.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["run", "url", "md5", "dest"])
            w.writerows(dl)
        if split:
            with open(OUTDIR / f"{study}.split.tsv", "w", newline="") as f:
                w = csv.writer(f, delimiter="\t", lineterminator="\n")
                w.writerow(["run", "src", "n_records", "out1", "out2"])
                w.writerows(split)
        with open(OUTDIR / f"{study}.samplesheet.csv", "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
            w.writerows(sheet)
        print(f"{study}: {len(rs)} runs, {len(dl)} files, {len(layouts)} samples, "
              f"{total / 1e9:.1f} Gbases", file=sys.stderr)


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["--all"]:
        args = sorted({r["study_id"] for r in csv.DictReader(open(ROOT / "metadata" / "runs.tsv"), delimiter="\t")
                       if r["include"] == "TRUE"})
    main(args)
