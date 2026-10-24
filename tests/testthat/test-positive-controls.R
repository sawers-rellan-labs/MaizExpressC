# Cold positive controls (HANDOFF M5/M8; annotation/cold_positive_controls.tsv): the eight maize DREB1/CBF
# genes are strongly induced by acute cold in leaves and shoots, not in roots or by chronic cool growth.

ctl <- tsv("annotation", "cold_positive_controls.tsv")
dreb <- ctl[role == "positive_control", gene_id]
meta <- pq("results", "meta_tissue_class.parquet")[filter_set == "all"]
cons <- tsv("metadata", "contrasts.tsv")

test_that("the control list holds 8 DREB1/CBF genes, 4 dehydrins and 3 ICE1 homologues", {
  expect_equal(ctl[, .N, by = role][order(role), N], c(8L, 3L, 4L))   # positive_control, reference, weak_control
  # seven have DDF2/DREB1F as Arabidopsis best hit; ereb23 has none in Ensembl (EnTAP: Setaria DREB1A) and sits
  # 7 kb from ereb36 in a tandem DREB1 cluster
  expect_equal(ctl[role == "positive_control" & ath_symbol == "ddf2", .N], 7L)
  expect_equal(ctl[role == "positive_control" & ath_symbol != "ddf2", symbol], "ereb23")
})

test_that("every DREB1/CBF control is strongly up in pooled leaf chilling", {
  m <- meta[gene_id %in% dreb & tissue_group == "Leaf" & contrast_class == "chilling"]
  expect_equal(nrow(m), 8)
  expect_true(all(m$k_studies >= 10))
  expect_true(all(m$pooled_log2fc > 3))
  expect_true(all(m$padj < 0.05))
})

test_that("most acute leaf and shoot contrasts induce at least half of the DREB1 controls", {
  de <- pq("results", "de.parquet")[gene_id %in% dreb & replicated == TRUE]
  acute <- cons[pool == "TRUE" & contrast_class %in% c("chilling", "freezing") &
                  tissue_group %in% c("Leaf", "Shoot", "Whole seedling") & as.numeric(cold_C) <= 10, contrast_id]
  up <- de[contrast_id %in% acute, .(up = sum(padj < 0.05 & log2fc > 0, na.rm = TRUE)), by = contrast_id]
  expect_gt(nrow(up), 60)
  expect_gt(mean(up$up >= 4), 0.7)
})

test_that("roots and chronic cool growth do not induce the DREB1 controls", {
  r <- meta[gene_id %in% dreb & tissue_group == "Root" & contrast_class == "chilling"]
  expect_false(any(r$padj < 0.05 & r$pooled_log2fc > 1))
  de <- pq("results", "de.parquet")[gene_id %in% dreb & contrast_id %in% cons[contrast_class == "chronic", contrast_id]]
  expect_equal(de[padj < 0.05 & log2fc > 1, .N], 0)
})

test_that("Li2022c wild type, with its deposit labels swapped back, induces DREB1 like the mutant", {
  de <- pq("results", "de.parquet")[gene_id %in% dreb]
  for (cid in c("Li2022c_chilling_LH244_Leaf", "Li2022c_chilling_bzip68-1_Leaf"))
    expect_gt(de[contrast_id == cid & gene_id == "Zm00001eb318890", log2fc], 5, label = cid)   # ereb36
  smp <- tsv("metadata", "samples.tsv")
  expect_equal(smp[sample_id == "Li2022c_WT_25C_12h_R1", phase], "cold")                  # titled 25 C, cold in fact
  expect_match(smp[sample_id == "Li2022c_WT_25C_12h_R1", curation_note], "CONDITION SWAPPED")
})
