#!/usr/bin/env python3
"""Write AWS SRA download manifests for the studies fetched from NCBI's AWS mirror (M3).

ENA delivers 2-8 MB/s per job from Hazel; NCBI's open-data bucket on AWS (sra-pub-run-odp)
gave ~50 MB/s per connection (21 Oct 2026). For the large studies in SRA_ROUTE, runs are
downloaded as .sra files (download_sra.sh, xfer node) and converted to FASTQ on compute
nodes (convert_sra.sbatch). The FASTQ paths match the study's nf-core samplesheet, so
run_study.sh is unchanged.

  pipeline/manifests/<study>.sra.tsv  run, url, out1, out2, reads_per_file

reads_per_file is the expected record count in each FASTQ: ENA read_count for ordinary
paired runs; half of it for single-file paired runs (Zhou2024: SRA stores R1 and R2 as
separate spots, which fasterq-dump --split-files writes to _1 and _2).
Usage (project root): python3 pipeline/make_sra_manifests.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAN = ROOT / "pipeline" / "manifests"
SRA_ROUTE = ["Crow2020", "Frey2020", "Zhou2022", "Zhou2024"]
SPLIT_HALVES = {"Zhou2024"}

runs = {r["run_accession"]: r for r in csv.DictReader(open(ROOT / "metadata" / "runs.tsv"), delimiter="\t")
        if r["include"] == "TRUE"}
for study in SRA_ROUTE:
    sheet = list(csv.DictReader(open(MAN / f"{study}.samplesheet.csv")))
    out = []
    for row in sheet:
        run = Path(row["fastq_1"]).name.split("_")[0]
        r = runs[run]
        n = int(r["n_spots"])
        out.append([run, f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}", row["fastq_1"], row["fastq_2"],
                    n // 2 if study in SPLIT_HALVES else n])
    with open(MAN / f"{study}.sra.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "url", "out1", "out2", "reads_per_file"])
        w.writerows(out)
    print(f"{study}: {len(out)} runs")
