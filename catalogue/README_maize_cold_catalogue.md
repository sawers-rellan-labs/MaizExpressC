# Maize cold RNA-seq catalogue — README

Built 20 Oct 2026 (M1). Reviewed by Rubén on 20 Oct 2026: only the 50 core studies were kept, and 5 studies found through their datasets were added (Jiang2022, Zeng2024, Gao2024, Liao2026, Crow2020), for 55. Purpose: list published RNA-seq experiments in maize (Zea) that compare a low-temperature treatment with a control temperature, with the setup and data accessions needed to download and re-quantify the reads.

## Files
- `maize_cold_rnaseq_studies.tsv`: 55 studies, one row per paper.
- `maize_cold_not_core.tsv`: 12 studies left out at review (Rubén, 20 Oct 2026), each with the reason. Gao2024, first tentative, was re-added once its cold RNA-seq dataset was found.
  - 3 borderline: post-harvest storage of detached kernels; a cold soaking test before germination at 25 °C; a design with no control-temperature libraries.
  - 9 tentative: the design could not be confirmed from full text (abstract only, or a preprint or journal outside PubMed).
- `maize_cold_candidate_additions.tsv`: what the search of unmatched repository records found, and the suggestion for each. Rubén added the 5 papers and left out the datasets without a paper.
- `maize_cold_related.tsv`: 63 related items. These are microarray and small-RNA studies, papers that only reanalyse public cold data (useful pointers to datasets), and cold transcriptome papers whose design was ambiguous.
- `maize_cold_to_check.tsv`: 10 open leads. These are dataset accessions reused by other papers whose primary study still has to be identified, plus surveys and meta-analyses to cross-check coverage against.
- `build/`: how the catalogue was made; see `build/search_log.md`.

## Inclusion
New transcriptome sequencing of Zea tissue under low temperature and a control temperature or a pre-cold baseline. That covers chilling, freezing, cold time courses, recovery, low-temperature germination and field cold. It also covers cold crossed with genotype (tolerant vs sensitive), mutants or transgenics, another stress, or a priming treatment, provided temperature itself was varied. Multi-stress panels are included when cold has its own control. The full rules are in `build/agent_instructions.md`.

## Key columns
- `design_category`: CHILL, FREEZE, TIMECOURSE, RECOVERY, ACCLIM, GERM, FIELD, CxGENO, CxSTRESS, MUTANT, CxTREAT, MULTI (joined with +).
- `cold_treatments`, `control_temperature`, `timepoints_duration`: as reported (°C, h or d); NR when not reported.
- `accession`: accessions the screening agents read in the paper, or found by web search; `notes` says which.
- `accession_harvested`: accessions found by `metadata/build/harvest_accessions.py` in the PMC and Europe PMC full-text XML (which keeps the data-availability sections) and through NCBI PubMed links, with their source.
- `accession_status`:
  - `stated`: in the paper per the agent.
  - `found by harvester`: only found by the harvester.
  - `lead`: unverified.
  - `to resolve`: none found yet.
- `related_to`: papers that share a dataset (count libraries once).
- `fit` (all `core` in the main table) and `fit_note`.

## Counts (55 studies; accession status as in the catalogue, see metadata/accessions_resolved.tsv for the resolved state)
- 24 with an accession stated in the paper, 10 found by the harvester, 1 lead, and 15 still to resolve. Most of the unresolved are recent papers whose data statements the full-text services drop, or papers without PMC full text.
- Repositories: mostly NCBI (BioProject/SRA/GEO). CNCB GSA holds Tang2026 and He2023a; ENA holds Bozic2025.
- Designs: chilling 47, tolerant vs sensitive or other genotype contrasts 27, time courses 14, germination 12, multi-stress panels 9, mutants 10, recovery 6, priming or chemical × cold 5, freezing 1.
- Tissues: leaf 29, root 13, seed/germination 10, shoot (including mesocotyl and coleoptile) 9, whole seedling 8.
- Years: 2015–2026, 46 of them since 2020.
- 5 core studies were screened from the abstract alone (the design was clear there); their full text is read at accession resolution.

## Known gaps
- Keyword searches miss papers whose abstracts never say "cold" or "RNA-seq". Makarevitch2015, the B73/Mo17/Oh43 cold, heat, salt and UV panel, is one of them; it was found through its dataset. The dataset-to-paper step (`build/dataset_papers.py`) covers maize cold BioProjects and GEO series that NCBI links to a paper, but not datasets without a PubMed link. Those will show up as unmatched repository records at accession resolution.
- Shared datasets: Zhang2024/Dou2025, Zhao2024a/b and the B144 vs Q319 papers (Yu2021, Yu2023, Yu2024) probably reuse libraries. This gets checked at accession resolution by comparing run lists.

## Next steps (after Rubén's review)
1. Accession resolution, as in MaizExpressN: verify every accession (organism, strategy, runs) with `verify_accessions.py`, and resolve the remaining 15. Match the unlinked repository records (`metadata/build/repo_*_maize_cold.tsv`) and the leads in `maize_cold_to_check.tsv`. Checkpoint: `metadata/accessions_resolved.tsv`.
2. Run tables and curated sample sheets (M2).
