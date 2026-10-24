# MaizeExpressC

Look up a maize gene and see how its expression responds to cold (chilling, freezing, cold during
germination, growth in a cool environment and recovery) across published RNA-seq studies, study by study and
pooled across studies by tissue, contrast class, temperature, time, growth system, stage and genotype.
GEMMA Lab, NC State (Rubén Rellán-Álvarez).

Sibling of [MaizExpressN](https://github.com/sawers-rellan-labs/MaizExpressN) (nitrogen), whose pipeline,
analysis code, gene annotation and app this project reuses. The plan is in `HANDOFF.md`; working notes and
every decision taken are in `CLAUDE.md`.

Large files (Parquet tables, per-study gene tables, MultiQC reports, annotation sources) are stored with
[Git LFS](https://git-lfs.com). Install it before cloning (`brew install git-lfs`, then `git lfs install`), or clone
with `GIT_LFS_SKIP_SMUDGE=1` to get only the code and pointer files.

## In numbers

34 studies, 1,148 libraries in use after library QC, 238 contrasts (162 replicated and tested, 76 shown
only), 35,008 genes tested. Classes: chilling, freezing, recovery, germination and chronic (grown cool);
eight DREB1/CBF genes serve as positive controls (`annotation/cold_positive_controls.tsv`).

## Run the app

```bash
Rscript -e 'shiny::runApp("app", launch.browser = TRUE)'
```

Run from the project root (the renv library activates there). The app reads `app/data/`, built by
`app/prepare_data.R` and `app/prepare_partners.R`; set `MAIZEEXPRESSC_DATA` to use a data folder elsewhere.

## Deploy to Posit Connect Cloud (not deployed yet)

Same route as MaizExpressN (account `gemmalab`). From an R session started in the project root:

```r
rsconnect::connectCloudUser()                                   # once per computer: authorise in the browser
rsconnect::deployApp(appDir = "app", appName = "maizeexpressc")  # upload app/ (about 720 MB) and publish
```

The upload is about 716 MB (`app/data` 721 MB: DE 383, TPM 206, partner lists 66, atlas 16), under the
1 GiB bundle limit of Connect Cloud's free plan; publishing from GitHub has no bundle limit. For that route,
refresh `app/manifest.json` with `rsconnect::writeManifest(appDir = "app")` and commit it with the `app/` folder.
Re-run `deployApp()` after `Rscript app/prepare_data.R` and `Rscript app/prepare_partners.R` to update the app.
Do a public build (`MAIZEEXPRESSC_PUBLIC=1`) if any study is marked private.

## Layout

| Folder | What is in it |
|---|---|
| `catalogue/` | The study catalogue (`maize_cold_rnaseq_studies.tsv`), the source of truth for which studies are in scope |
| `metadata/` | Curated tables: `accessions_resolved.tsv`, `runs.tsv`, `samples.tsv`, `contrasts.tsv`, `needs_review.tsv`, `library_qc.tsv`, `library_exclusions.tsv`; builders and per-study sample sheets in `metadata/build/` |
| `pipeline/` | Download and quantification on NC State's Hazel cluster (Slurm, nf-core/rnaseq; ENA downloads, or AWS SRA with conversion for the largest studies); see `pipeline/README.md` |
| `results/` | Per-study gene tables copied from Hazel (`results/<study>/`), `de.parquet` and `de_summary.tsv` (M4), `meta_*.parquet` (M6), `coexpr_strongest_pairs.tsv` |
| `annotation/` | Gene annotation (M5, as in MaizExpressN): sources in `annotation/source/`, gene, alias, GO, MaizeMine, literature and atlas tables, and the cold positive controls |
| `analysis/` | `de_per_study.R` (DESeq2 per study, with replicate blocking and the paired diversity-panel fit) and `pool_meta.R` (pooling across studies) |
| `app/` | The Shiny app (`app.R`, `R/`, `www/`, `methods.md`), `prepare_data.R` and `prepare_partners.R`. Pages: Gene lookup, Gene list, Co-expression, Rank genes, Studies, Methods |
| `tests/` | testthat suite: `Rscript tests/testthat.R` |

## Rebuild everything after a change

From the project root, in this order (each step reads the outputs of the one before):

```bash
python3 metadata/build/build_runs.py
python3 metadata/build/merge_samples.py --write
python3 metadata/build/build_contrasts.py
python3 metadata/build/validate_contrasts.py
python3 metadata/build/library_exclusions.py
Rscript analysis/de_per_study.R            # ~30 min on 6 cores (DE_WORKERS)
Rscript analysis/pool_meta.R
Rscript annotation/build_parquet.R         # after the annotation fetchers, if the sources change
Rscript annotation/build_cold_controls.R
Rscript app/prepare_data.R
Rscript app/prepare_partners.R             # co-expression partner lists, ~17 min
Rscript tests/testthat.R
```

## Add a study

1. **Catalogue.** Add a row to `catalogue/maize_cold_rnaseq_studies.tsv` with a new `study_id`
   (FirstAuthorYear). For unpublished lab data add a column `private` and set it to `TRUE`, then build the
   pooled tables and the app data with `MAIZEEXPRESSC_PUBLIC=1` for a public deployment.
2. **Accessions and runs.** Record the accession in `metadata/build/resolve_accessions.py` (→
   `metadata/accessions_resolved.tsv`), then run `build_runs.py`; run-level exclusions go in its rules with a reason.
3. **Samples.** Write `metadata/build/samples/<study_id>.tsv` (and a `.review.tsv`) following
   `metadata/build/curation_spec.md`: genotype, cold tolerance, tissue, temperature (day and night), control
   temperature, phase (control / cold / recovery / chronic), time since the cold started and recovery time.
   `merge_samples.py --write` validates it.
4. **Contrasts.** `build_contrasts.py` pairs each cold group with its control (time-matched, else the single
   control, else time 0) and assigns the class; add a study rule only for special designs. Check with
   `validate_contrasts.py`.
5. **Quantify on Hazel.** Make the manifests, sync `pipeline/` to Hazel, download and run the head job
   (details in `pipeline/README.md`). Copy `salmon.merged.gene_counts.tsv`, `gene_tpm.tsv`, `gene_lengths.tsv`
   to `results/<study_id>/`.
6. **QC.** Run `pipeline/library_qc.py` on Hazel into `metadata/library_qc.tsv`, then `library_exclusions.py`.
   Check the DREB1/CBF controls respond where they should (leaf or shoot cold at 10 °C or below).
7. **Analyse and publish.** `de_per_study.R`, `pool_meta.R`, `app/prepare_data.R`, `app/prepare_partners.R`, then the tests.

## Conventions

- Gene key: B73 RefGen_v5 gene IDs (`Zm00001eb…`); v4, v3, symbols and Arabidopsis names are search aliases.
- Effect sign: log2FC = log2(cold / control), so positive means higher under cold; recovery contrasts are
  log2(recovered / control). The app's sign toggle flips it.
- Counts are never pooled across studies; only log2FC with its SE. Chilling and freezing pool together as
  "acute cold"; germination, chronic and recovery designs stay apart.
- Library QC is strict (`metadata/build/library_exclusions.py`); unreplicated contrasts are shown but never pooled.

## License

The code is released under the [MIT License](LICENSE) (© 2026 Rubén Rellán-Álvarez and the GEMMA Lab).
The expression results are derived from published studies whose raw reads are in public repositories;
please cite the original studies (listed with their DOIs in `catalogue/maize_cold_rnaseq_studies.tsv`) when using them.
