# metadata/build — accession resolution (M1)

Provenance for `metadata/accessions_resolved.tsv`. Python 3 standard library only; responses are cached in `cache/` (git-ignored), so reruns are offline and fast. No sequence data is downloaded. Run each script from the project root, in this order.

| Script | Output | What it does |
|---|---|---|
| `search_repositories.py` | `repo_bioproject_maize_cold.tsv`, `repo_geo_maize_cold.tsv`, `repo_sra_maize_cold.tsv` | All maize BioProjects, GEO series and SRA studies that mention cold, for matching. |
| `harvest_accessions.py` | `accession_evidence.tsv` | For each catalogue study: NCBI elink from PMID to BioProject, SRA and GEO; Europe PMC annotations and database links; accession patterns in the PMC and Europe PMC full-text XML. |
| `build_candidates.py` | `accession_candidates.tsv` | Merges the catalogue's accessions, the harvester's, and those found by hand (`EXTRA`, with their source). |
| `verify_accessions.py` | `accession_records.tsv`, `ena_runs/<acc>.tsv` | Fetches each record: BioProject / GEO summary, plus the ENA run report (organism, strategy, layout, instrument, read counts, FASTQ URLs and MD5s). |
| `resolve_accessions.py` | `../accessions_resolved.tsv` | One row per study and accession, with status, role, evidence and maize RNA-seq run count; the decisions are in its `D` table. |
| `unmatched_records.py` | `uncatalogued_repository_records.tsv` | Maize cold records that match no catalogue study (`likely_cold_rnaseq` marks the plausible ones). For Rubén to triage; nothing is added to the catalogue. |

Steps done by hand (recorded in the `evidence` and `notes` columns):
- Data statements read on PMC web pages. These keep back matter that the PMC XML drops; this is how Zhao2021's SRA runs and the not-deposited statements were found.
- Publisher pages for studies without PMC text, read in a browser: OUP (Zhou2024), Wiley (Waters2017, Meng2025) and Springer (Han2026). ScienceDirect (Friero2023) did not load.
- CNCB GSA records (CRA007319, CRA027029, CRA041093, CRA041107) checked in a browser, because the GSA site blocks scripted requests.
- Line-name searches in NCBI (e.g. B144/Q319, He344, 18-599) and submitter/design matches against the repository tables. One match (Han2026 → PRJNA1314463) was refuted by the paper's own statement.

## M2: run tables and sample sheets

| Script / file | Output | What it does |
|---|---|---|
| `build_runs.py` | `metadata/runs.tsv` | One row per run of every resolved dataset: ENA run report, FASTQ URLs and MD5s, and the ENA sample record (submitter alias, title, attributes) as JSON. Include/exclude rules are written per study in the script, each with a reason. Shared libraries get one owning `study_id` plus `also_used_by`. |
| `gsa_runs.tsv`, `gsa_md5_*.txt` | (input) | CNCB GSA run and sample metadata captured by hand from the GSA web pages (the site blocks scripts), and the MD5 lists from download.cncb.ac.cn. |
| `curation_spec.md` | (instructions) | Field definitions and vocabularies for the per-study sample sheets. |
| `samples/<study_id>.tsv`, `.review.tsv` | (curated) | One sample sheet and one list of open questions per study. |
| `review_global.tsv` | (curated) | Project-wide open questions. |
| `merge_samples.py` | `metadata/samples.tsv`, `metadata/needs_review.tsv` | Validates the sheets against runs.tsv (each included run exactly once, vocabularies, numbers) and merges them with `--write`. |

## M4 (first step): contrasts

| Script | Output | What it does |
|---|---|---|
| `build_contrasts.py` | `metadata/contrasts.tsv` | Groups included samples into condition groups and builds contrasts per study. Generic rule: within strata of held-constant factors (genotype, tissue, other factor, time), every lower-N group vs every higher-N group. Per-study rules (class, strata, exclusions, notes, and the special designs: Jia2023, Braynen2025, Wang2019, Ravazzolo2020 form, Li2026a split-root, Meng2020, Wang2020a) are in the `rules()` table. |
| `validate_contrasts.py` | (report) | Checks sample ids, disjoint groups, n counts, sign convention, held-constant factors, duplicate pooled contrasts, and lists included samples no contrast uses. |

Columns beyond the HANDOFF schema: `numerator_samples` / `denominator_samples` (sample ids for DESeq2), `genotype`, `tissue_group`, `tissue_detail`, `other_factor_level` (for filters and pooling), `tags`, `pool` (FALSE = display only) and `batch_term` (a blocking factor to add to the study's DESeq2 design).

Tags: `unreplicated`, `single_library_arm`, `not_time_matched` (the control is a baseline at another time), `batch_confounded`, `batch_blocked`, `local_supply`, `systemic_signal`, `split_vs_homogeneous`, `merged_replicates`, `low_depth`, `unexplained_sample_sets`, `deferred_processing`, `label_check`.

`time_h` in contrasts is hours since that contrast's own N change (for a resupply after starvation, hours since resupply); in `samples.tsv` it is on the experiment clock.

## Status vocabulary in `accessions_resolved.tsv`

- `resolved`: the paper states the accession, or the repository record links the PMID, or (`evidence_type = repository_match`) the record's title, authors or design match the paper with no conflict. Run contents were checked against the design.
- `ambiguous`: there is a candidate, but the link to the paper is unconfirmed.
- `unavailable`: the authors say the data aren't deposited (not shared, on request, or in the supplement only).
- `not_found`: the paper gives no accession (or couldn't be read) and no repository record matched.

`role` says how a row's runs are used: `primary`; `shared` (the same runs are claimed by another study, count once); `partial` (only some runs are new); `duplicate` (a re-upload, don't use); `exclude` (out of scope: sorghum, small RNA, ChIP); `review` (in scope is unclear); `umbrella` (a GEO SuperSeries); `candidate` (goes with `ambiguous`).
