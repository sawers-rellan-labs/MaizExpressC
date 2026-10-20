# Search log (20 Oct 2026)

Provenance for `catalogue/maize_cold_rnaseq_studies.tsv`. Searches and screening were run by Claude Code for Rubén Rellán-Álvarez, following the MaizExpressN catalogue workflow.

## 1. PubMed (`pubmed_search.py`, NCBI E-utilities)

Ten queries; hits and new PMIDs per query are in `pubmed_queries.tsv`.

| # | Focus | Hits | New |
|---|---|---|---|
| 1 | maize/Zea/corn × cold/chilling/low temperature/freezing/frost × RNA-seq/transcriptome/DEGs | 325 | 325 |
| 2 | MeSH: Zea mays × Cold Temperature/Cold-Shock Response/Freezing × expression profiling/transcriptome/RNA-Seq | 155 | 83 |
| 3 | cold stress / chilling tolerance / cold acclimation / freezing tolerance phrases × transcriptome | 260 | 61 |
| 4 | cold × DEGs, WGCNA, co-expression, lncRNA, circRNA, single-cell, alternative splicing, eQTL | 119 | 20 |
| 5 | cold × germination / emergence / seed vigour × transcriptome | 103 | 20 |
| 6 | abiotic stress panels × RNA-seq × cold/temperature | 151 | 38 |
| 7 | suboptimal / low night temperature, early sowing, cool temperature × expression | 99 | 33 |
| 8 | cold × ICE1/CBF/DREB/cold-responsive genes × transcriptome | 20 | 0 |
| 9 | teosinte, mexicana, highland, temperate/tropical × cold × expression (restricted to Zea) | 479 | 127 |
| 10 | cold × microRNA/small RNA/degradome/microarray (for the related table) | 65 | 11 |

Union: 718 PMIDs (`pubmed_candidates.tsv`). A lenient keyword pre-screen set aside 356 records whose title and abstract lack a Zea term, any temperature term, or any expression term; they are listed with the reason in `pubmed_prescreen_excluded.tsv`. Records without an abstract went to screening. The other 362 were split into 11 batches (`batches/`).

## 2. Screening and extraction

Each batch was screened by a subagent with the rules in `agent_instructions.md`. Abstracts were triaged, then included papers were read in full text from PubMed Central where available; outputs are in `agent_out/pubmed_batch_*.tsv`. The agents' RELATED items are collected in `agent_out/related_reports.tsv`, and the leads they raised in `agent_out/to_check.tsv`.

## 3. Beyond PubMed

- **Dataset-to-paper links (`dataset_papers.py`):**
  - Inputs: every maize BioProject, GEO series and SRA study that mentions cold (`metadata/build/repo_*_maize_cold.tsv`, from `search_repositories.py`: 122 BioProjects, 37 GEO series and 84 SRA projects), plus the dataset leads from the agents.
  - These were linked to PubMed with NCBI elink. Of the 66 papers found, 47 had not been screened yet; with Makarevitch2015 (PRJNA244661, not linked to its paper in NCBI) added by hand, they formed `batches/batch_dataset.tsv`. Its screening agent included Makarevitch2015.
- **Web search** (preprint servers, journals PubMed indexes poorly, gene-function papers, expression atlases): `agent_out/web_literature.tsv`. It added Mao2017 (Crop Journal), Wang2025a (Agronomy), and the preprints Cackett2023 and James2026, all tentative. Its queries are listed in the agent's report (summarised in the session notes).
- **Manual addition:** the 11-angiosperm cold time course descriptor (Wang2023, Scientific Data, PRJNA767196). NCBI shows 12 Zea mays runs in it (`agent_out/manual_additions.tsv`).

## 4. Accessions

`metadata/build/harvest_accessions.py` (pointed at the cold catalogue) collected candidate accessions per study: NCBI elink from PMID to BioProject, SRA and GEO; Europe PMC annotations; and accession patterns in the PMC and Europe PMC full-text XML. The XML keeps the data-availability sections that the agents' full-text tool drops. Results are in `metadata/build/accession_evidence.tsv` and the catalogue's `accession_harvested` column.

## 5. Merge

`build_catalogue.py` merges the agent tables and applies manual curation (the `OV` table): study-ID clashes, shared datasets, and borderline or tentative fits. It writes the catalogue, related and to-check tables.

## Known limits

- The PubMed full-text tool used by the agents drops data-availability sections; the harvester recovers many of them, but 15 core studies still have no accession.
- Publisher pages (Elsevier, Springer, Wiley, MDPI) often blocked scripted fetches, so studies without PMC full text are abstract-level.
- Keyword recall: papers that never mention cold or RNA-seq in the abstract are found only through their datasets (section 3).
