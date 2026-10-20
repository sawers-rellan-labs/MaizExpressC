# MaizeExpressC

Look up a maize gene and see how its expression responds to cold (chilling and freezing) across published
RNA-seq studies, study by study and pooled across studies by tissue, contrast class, growth system, stage and
genotype. GEMMA Lab, NC State (Rubén Rellán-Álvarez).

Sibling of [MaizExpressN](https://github.com/sawers-rellan-labs/MaizExpressN) (nitrogen), whose pipeline,
analysis code, gene annotation and app this project reuses. The plan is in `HANDOFF.md`; working notes and
decisions are in `CLAUDE.md`.

## Status

M0 (scaffold) done. M1: `catalogue/maize_cold_rnaseq_studies.tsv` lists 50 maize cold RNA-seq studies (reviewed 20 Oct 2026; see `catalogue/README_maize_cold_catalogue.md`). Next: accession resolution.

## Layout

| Folder | What is in it |
|---|---|
| `catalogue/` | The study catalogue and how it was built (`catalogue/build/`) |
| `metadata/` | Accessions, runs, samples and contrasts once curated; builders in `metadata/build/` |
| `pipeline/` | Download and quantification on NC State's Hazel cluster (Slurm, nf-core/rnaseq); see `pipeline/README.md` |
| `analysis/` | `de_per_study.R` (DESeq2 per study) and `pool_meta.R` (pooling across studies) |
| `annotation/` | B73 RefGen_v5 gene annotation (same as MaizExpressN) |
| `app/` | The Shiny app (template from MaizExpressN until M7) |
| `tests/` | testthat suite: `Rscript tests/testthat.R`; N tests to adapt are in `tests/templates_from_N/` |

## Conventions

- Gene key: B73 RefGen_v5 gene IDs (`Zm00001eb…`); v4, v3, symbols and Arabidopsis names are search aliases.
- Effect sign (proposed): log2FC = log2(cold / control temperature), positive = higher under cold.
- Counts are never pooled across studies; only log2FC with its SE.
