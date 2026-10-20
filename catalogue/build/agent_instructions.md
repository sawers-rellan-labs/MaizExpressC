CONTEXT
You are helping a maize geneticist (NC State, GEMMA lab) build a catalogue of published RNA-seq experiments in maize grown under cold (chilling or freezing) versus a control temperature. The catalogue will drive bulk download and uniform re-quantification of raw reads, so the EXPERIMENTAL DESIGN (temperatures, durations, genotypes, tissues, replicates) and DATA ACCESSIONS matter most. Accuracy over speed: never guess.

INPUT
A TSV batch of PubMed records (pmid, year, journal, title, doi, pmcid, publication_types, abstract). Read it with the Read tool.

STEP 1 - TRIAGE (every record in the batch)
INCLUDE a paper only if BOTH:
 (a) it generated NEW transcriptome sequencing data: bulk mRNA-seq, total/strand-specific RNA-seq (incl. lncRNA/circRNA/alternative-splicing studies), 3' tag-seq/QuantSeq, cell-type/LCM RNA-seq, single-cell/single-nucleus or spatial RNA-seq, or Iso-Seq/full-length combined with expression; AND
 (b) maize (Zea mays; teosinte/other Zea is OK but flag it) was sampled under a low-temperature treatment AND a control temperature (or a time-0/pre-cold baseline): chilling (above 0 °C), freezing (0 °C or below, with or without acclimation), cold time courses, recovery after cold, low-temperature germination/emergence vs optimal temperature, field cold (early sowing or a natural cold episode) with a warmer control, or cold crossed with genotypes (tolerant vs sensitive), mutants/transgenics, another stress, or a chemical/priming treatment PROVIDED temperature itself was varied. Multi-stress panels (e.g. cold, heat, drought, salt each vs one control) are INCLUDED when cold has a control.
Put in RELATED (not the main table): microarray-only; small RNA/miRNA/degradome-only; ribosome-profiling-only; papers that only REANALYSE public RNA-seq data (record which datasets they used, if stated); cold transcriptome papers where the abstract is ambiguous and you could not confirm the design from full text.
EXCLUDE: qPCR-only gene-function papers (overexpression/knockout of one gene with qPCR or no transcriptomics), proteomics/metabolomics-only, reviews, heat-only or other stresses without a cold treatment, non-Zea species (also when maize is only a gene donor expressed in Arabidopsis/rice/tobacco), seed storage/germination studies without a temperature contrast.
If an abstract is ambiguous (e.g. "transcriptome" without saying RNA-seq vs microarray, or cold not clearly a treatment), check the full text before deciding.

STEP 2 - EXTRACT (each included study)
Load the PubMed tools first with ToolSearch: query "select:mcp__98767f22-c23e-4bdc-88f3-773f464b12fe__get_full_text_article,mcp__98767f22-c23e-4bdc-88f3-773f464b12fe__convert_article_ids,mcp__98767f22-c23e-4bdc-88f3-773f464b12fe__get_article_metadata". Read the full text (one article per call) when the batch gives a PMCID: Methods, figure legends, Data Availability / accession numbers. If there is no PMC full text, use the abstract plus WebSearch (e.g. the paper title + "BioProject" / "PRJNA" / "GSE" / "SRA" / "CRA") to find the accession; at most ~3 WebFetch calls in total, only for a publisher page when essential. If a fetch is refused or rate-limited, do not retry it and do not try other routes (no curl, Python, mirrors or caches): mark the field 'not found'. Never use Bash for web access.
NEVER invent or infer accessions, numbers, units, temperatures or genotypes. Use 'NR' when a detail is not reported and 'not found' for an accession you could not locate. Distinguish 'not stated' (no data statement) from 'on request'.

STEP 3 - OUTPUT
Write a TSV (tab-separated, UTF-8, header row, no tabs or newlines inside fields; use '; ' inside fields) with the Write tool to the path you were given. Do not modify any other files. Columns, in this order:
study_id	pmid	pmcid	doi_url	first_author	year	journal	title	design_category	genotypes	n_genotypes	tissue	stage_age	growth_system	cold_treatments	control_temperature	other_factors	timepoints_duration	bio_reps	n_libraries	platform_layout	accession	repository	other_omics	key_finding	source_of_details	notes
Field guidance:
- study_id: FirstAuthorSurnameYear (e.g. Li2020); add a/b if needed.
- doi_url: https://doi.org/<doi>
- design_category: one or more of CHILL (cold above 0 °C vs control), FREEZE (0 °C or below), TIMECOURSE (several cold durations), RECOVERY (back to control temperature after cold), ACCLIM (cold acclimation before freezing), GERM (low-temperature germination/emergence), FIELD (field cold or early sowing), CxGENO (contrasting genotypes/hybrids/panels), CxSTRESS (with another stress), MUTANT (mutant/transgenic x cold), CxTREAT (chemical/priming/hormone x cold), MULTI (cold is one of several stresses vs a shared control); join with '+'.
- genotypes: names (e.g. 'B73; Mo17'); note claimed cold tolerance (tolerant/sensitive) and flag teosinte or highland material.
- tissue: e.g. 'seedling leaf', 'shoot apex', 'root tip', 'germinating seed/embryo', 'mesocotyl'.
- stage_age: e.g. 'V3', '14 d after sowing', 'three-leaf stage'.
- growth_system: growth chamber / greenhouse / field (+ location, years); soil pots / hydroponics / paper rolls; light and photoperiod during cold if stated.
- cold_treatments: temperature(s) in °C with day/night if given, e.g. '4 °C', '10/8 °C day/night', '-4 °C for 2 h after 7 d at 4 °C (acclimated)'.
- control_temperature: e.g. '25 °C', '28/22 °C', 'time 0 before cold'.
- other_factors: e.g. 'tolerant vs sensitive inbreds', 'drought', 'ABA priming', 'mutant zmice1', or 'none'.
- timepoints_duration: cold duration and sampling times, recovery times.
- bio_reps: biological replicates per condition (as stated). n_libraries: total RNA-seq libraries if stated or computable.
- platform_layout: e.g. 'Illumina HiSeq 2500, PE150', 'BGISEQ-500 SE50'.
- accession: ALL accessions (BioProject PRJNA/PRJEB/PRJDB/PRJCA, SRA SRP/SRR, GEO GSE, CNCB GSA CRA, ArrayExpress E-MTAB), '; '-separated; or 'on request' / 'not stated' / 'not found'.
- repository: NCBI SRA / NCBI GEO / CNCB GSA / ENA / DDBJ / Figshare / other.
- other_omics: metabolomics, proteomics, small RNA, physiology, GWAS... or 'none'.
- key_finding: <=20 words.
- source_of_details: 'PMC full text', 'abstract + web search', 'publisher page', etc.
- notes: caveats (e.g. 'pooled samples, no replicates', 'same dataset as <study_id>', 'teosinte').

FINAL REPLY (compact - the TSV file is the main deliverable):
 1. Path of the TSV and number of included studies.
 2. RELATED: one line each - 'PMID | FirstAuthorYear | data type (microarray/sRNA/degradome/reanalysis/ambiguous) | accession(s) if stated | <=10-word description'.
 3. EXCLUDED: comma-separated 'PMID (<=4-word reason)'.
 4. Problems (tools that failed, papers you could not read).
If writing the file fails, paste the TSV in the reply instead.
