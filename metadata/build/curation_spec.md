# Sample-sheet curation spec (M2, cold)

Goal: for each study, map every run in `metadata/runs.tsv` (rows with that `study_id`) to a biological sample with its experimental design. One study's output goes in `metadata/build/samples/<study_id>.tsv`; open questions go in `metadata/build/samples/<study_id>.review.tsv`. Write only those two files for your studies.

Accuracy over speed. Never guess a design detail: if the run metadata and the paper don't settle it, use `NR` (not reported) or flag it in the review file. Don't download sequence data.

## Sources, in order of trust
1. Run and sample metadata in `metadata/runs.tsv`: `sample_title`, and in `biosample_attributes` (JSON) the submitter `alias`, `title` and attributes (treatment, tissue, cultivar, time, replicate...). These say which run is which sample.
2. The paper. For PMC papers use the PubMed full-text tool, or the PMC web page (`https://pmc.ncbi.nlm.nih.gov/articles/<PMCID>/`), which also keeps back matter. Use Methods, figure legends and tables. GEO sample pages (`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM...&targ=self&form=text`) often add characteristics.
3. The catalogue row (`catalogue/maize_cold_rnaseq_studies.tsv`): a summary written by an earlier pass. Helpful but may be wrong; the paper wins.
4. `metadata/accessions_resolved.tsv` notes for that study (shared runs, partial projects, exclusions).

Be polite to servers: at most one request per second.

## Output: `<study_id>.tsv` (tab-separated, UTF-8, header row, no tabs or newlines in fields)

One row per biological sample (library). Several runs of the same library (lanes, re-sequencing) belong to one sample. Include every run of the study that has `include = TRUE` in runs.tsv, exactly once. Runs with `include = FALSE` are left out.

| Column | Content |
|---|---|
| sample_id | `<study_id>_<label>`, e.g. `Li2020a_4C_R1`. Unique, no spaces. |
| study_id | as in the catalogue |
| run_accessions | `;`-joined run accessions for this sample |
| genotype | as named in the paper, e.g. `B73`, `Mo17`, `zmice1-1`, `B73xMo17` |
| genotype_class | one of `B73`, `other inbred`, `hybrid`, `mutant`, `transgenic`, `teosinte`, `NR` (mutant = induced or natural mutant line incl. CRISPR knockouts; transgenic = overexpression / RNAi; introgression or near-isogenic lines = `other inbred`) |
| background | genetic background for mutants/transgenics (e.g. `B73`, `LH244`, `B104`); the line itself for inbreds/hybrids |
| cold_tolerance | the paper's own label for this genotype: `tolerant`, `sensitive`, `intermediate`, or `NR` |
| tissue_group | one of `Root`, `Leaf`, `Shoot`, `Whole seedling`, `Seed/germination`, `Reproductive` (Shoot = stem, shoot apex, mesocotyl, coleoptile or whole shoot; Seed/germination = imbibed or germinating seeds, embryos, radicles) |
| tissue_detail | as specific as the paper allows, e.g. `second leaf`, `primary root tip (0-1 cm)`, `root hairs`, `mesocotyl` |
| po_term | Plant Ontology ID if obvious (`PO:0025034` leaf, `PO:0009005` root, `PO:0020127` primary root, `PO:0009009` embryo, `PO:0000015` root hair cell); else empty |
| stage | as written, e.g. `V3`, `three-leaf stage`, `14 d after sowing`, `48 h imbibition` |
| stage_class | `seed/germination`, `seedling` (up to ~V4 or <= 21 d), `vegetative`, `reproductive`, or `NR` |
| age_days | plant age at sampling in days if stated or computable, else empty |
| growth_system | free text, e.g. `soil pots, growth chamber 14/10 h light`, `paper rolls, incubator`, `field, Jilin 2020` |
| temp_label_paper | label used in the paper/metadata, e.g. `CK`, `LT`, `cold`, `4C`, `T1`, `HF24h` |
| temperature_C | number only: the temperature the sample was at during its condition (day temperature when day/night differ); controls get the control temperature |
| temperature_night_C | night temperature if day/night differ, else empty |
| control_temperature_C | the study's control temperature (day value) for this sample's comparison |
| temp_rank_in_study | integer, 1 = coldest condition within the study (ties share a rank) |
| phase | `control` (never cold), `cold` (in cold at sampling, acute treatment), `recovery` (back at control temperature after cold), `chronic` (grown or germinated at that temperature from sowing / imbibition, no shift) |
| time_h | hours since cold onset, for every cold or recovery sample, and for controls of time-resolved designs (time-matched; 0 for time-zero samples). Empty for chronic designs and for single-point controls. |
| recovery_h | hours since return to the control temperature (recovery samples only), else empty |
| light_condition | free text if the paper says, e.g. `light`, `dark`, `16/8 h`, `continuous light during cold`; else empty |
| other_factor | non-temperature factor varied in the design, e.g. `hormone treatment`, `priming`, `humidity`, `tissue zone`, `acclimation`; empty if none |
| other_factor_level | e.g. `EBR 0.1 µM`, `methylglyoxal`, `high soil moisture`, `acclimated 7 d at 4 °C`, `control` |
| replicate | biological replicate number within its group (1, 2, 3...) |
| pooled_plants | plants pooled per sample if stated (a single number; ranges go in curation_note), else empty |
| include | `TRUE`, or `FALSE` with the reason in curation_note (e.g. failed library, or a sample outside any cold contrast with no control) |
| curation_note | short provenance: how the run was mapped (e.g. `alias LT6h-L2`), plus anything unusual (batch differences, day/night details, mismatches with the paper) |

Notes on the design fields:
- Every sample must be placeable in a cold contrast later: cold vs its control (same genotype, tissue, other factor and, for time courses, time point when the design has time-matched controls), recovery vs cold or vs control, chronic cool vs chronic warm.
- Freezing is a temperature at or below 0 °C; acclimation before freezing goes in other_factor (`acclimation`).
- Low-temperature germination (seeds germinated at, e.g., 10 °C from imbibition vs 25 °C): phase `chronic`, stage_class `seed/germination`; give the imbibition duration in stage.
- For mutant or transgenic x cold designs, both genotypes get their own samples; temperatures stay in the temperature fields.
- Shared runs: if a run's `also_used_by` lists another study, curate it under its owning `study_id` only.
- Samples that are not part of any cold contrast (e.g. a hormone-only treatment at the control temperature) are still listed with `include = TRUE` and explained in curation_note; contrasts decide later. Use `include = FALSE` only for samples that clearly can't be used.

## Output: `<study_id>.review.tsv`

Columns: `study_id`, `item`, `question`, `proposed_answer`, `evidence`. One row per open issue: an ambiguous run-to-sample mapping, a mismatch between paper and repository (replicate counts, labels), missing temperatures or durations, unexpected extra runs, batch confounding, dataset ownership doubts, anything a human must decide. Write the header even if there are no issues.

## Finish
Check your file: every included run of the study appears exactly once; columns are exactly as listed, in that order; controlled vocabularies are respected. Then reply with a 3-5 line summary per study (samples, design, open issues).
