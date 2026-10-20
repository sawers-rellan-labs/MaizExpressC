# Handoff: MaizeExpressC

Plan for a cold-stress counterpart of MaizExpressN: look up a maize gene and see how its expression responds to low temperature across published RNA-seq studies, study by study and pooled across studies by tissue, contrast class, growth system, stage and genotype. Owner: Rubén Rellán-Álvarez (GEMMA Lab, NC State). Started 20 Oct 2026 from the MaizExpressN infrastructure (see `CLAUDE.md` for what was copied).

## 1. Scope

In scope: new bulk RNA-seq (also lncRNA/total, 3′ tag-seq, LCM, single-cell, Iso-Seq with expression) of Zea tissue under at least two temperature conditions where one is low temperature:

- chilling (above 0 °C, typically 2–15 °C) vs a control temperature, steady or as a time course;
- freezing (0 °C or below), with or without cold acclimation;
- recovery after cold (return to the control temperature);
- cold at germination or emergence (low-temperature germination tests);
- field cold (early sowing, natural low temperatures) when a warmer control is sampled;
- cold crossed with genotype (tolerant vs sensitive lines, mutants, transgenics), with another stress, or with a treatment (e.g. a priming chemical) provided temperature itself is varied.

Related, not in the main table: microarray-only, small RNA/degradome-only, studies that only reanalyse public data, and multi-stress panels where cold is one of several treatments (these go in the main table when cold has its own control).

Out: heat-only, drought-only or other stresses without a cold treatment; qPCR-only; proteomics/metabolomics-only; non-Zea species.

## 2. Milestones (⏸ = review checkpoint: stop for Rubén)

| | Milestone | Reuses from MaizExpressN | New for cold |
|---|---|---|---|
| M0 | Repo scaffold | layout, pipeline, accession tools, annotation, renv | this plan, cold notes |
| M1 ⏸ | Catalogue: literature and repository search, then accession resolution | `catalogue/build` screening workflow, `metadata/build` NCBI/ENA tools | cold search terms and inclusion rules; checkpoints: catalogue, then `accessions_resolved.tsv` |
| M2 ⏸ | Run tables and curated sample sheets | `build_runs.py`, `merge_samples.py` | curation fields: temperature (°C), control temperature, duration, phase (steady / shock / recovery / acclimation), time_h, light/dark; checkpoint: `needs_review.tsv` |
| M3 | Download and quantification on Hazel | downloader, nf-core/rnaseq head job, shared B73 v5 salmon index | Rubén confirms storage paths first |
| M4 ⏸ | Contrasts and DESeq2 per study | `de_per_study.R`, strict library QC | cold classes and sign; checkpoint: `contrasts.tsv` |
| M5 | Gene annotation | `annotation/` as is | cold positive controls (e.g. CBF/DREB1, COR, ICE1 homologues) |
| M6 | Pooling across studies | `pool_meta.R` (DerSimonian–Laird) | class groups for cold |
| M7 ⏸ | Shiny app (first slice for review) | `app/` | labels, filters (temperature, duration, tolerance) |
| M8 | Tests, repository, deployment | testthat suite, Posit Connect Cloud route | GitHub repo once Rubén chooses visibility |

## 3. Proposed contrast classes (to settle at M4)

- `chilling`: cold above 0 °C vs control temperature.
- `freezing`: 0 °C or below vs control (acclimated or not, recorded).
- `recovery`: after return to control temperature vs still in cold (or vs control).
- `germination`: low-temperature germination or emergence vs optimal temperature.
- `field_cold`: field cold episode or early sowing vs warmer control.

Sign: log2FC = log2(cold / control), positive = higher under cold. Genotype tolerance (tolerant vs sensitive) is a filter and a stratum, not a class.

## 4. Layout

Same as MaizExpressN: `catalogue/`, `metadata/`, `pipeline/`, `analysis/`, `annotation/`, `results/`, `app/`, `tests/`, plus `private/` (gitignored) for anything that must not be public.
