# Contrast construction and the sign convention (HANDOFF M4): log2FC = log2(cold / control); recovery
# contrasts put the recovered group in the numerator.

cons <- tsv("metadata", "contrasts.tsv")
smp <- tsv("metadata", "samples.tsv")
setkey(smp, sample_id)

test_that("metadata/build/validate_contrasts.py finds no problem", {
  py <- Sys.which("python3"); skip_if(!nzchar(py), "python3 not found")
  out <- system2(py, file.path(ROOT, "metadata", "build", "validate_contrasts.py"), stdout = TRUE, stderr = TRUE)
  expect_match(out[1], "; 0 problems$")
})

test_that("every contrast refers to included samples of one tissue group, with disjoint sides", {
  for (i in seq_len(nrow(cons))) {
    c1 <- cons[i]; num <- split_ids(c1$numerator_samples); den <- split_ids(c1$denominator_samples)
    expect_true(length(num) > 0 && length(den) > 0, info = c1$contrast_id)
    expect_true(all(c(num, den) %in% smp$sample_id), info = c1$contrast_id)
    expect_true(all(smp[c(num, den), include] == "TRUE"), info = c1$contrast_id)
    expect_length(intersect(num, den), 0)
    expect_equal(unique(smp[c(num, den), tissue_group]), c1$tissue_group, info = c1$contrast_id)
  }
})

test_that("the colder side is the numerator; recovery contrasts compare recovered plants with controls", {
  for (i in seq_len(nrow(cons))) {
    c1 <- cons[i]; num <- smp[split_ids(c1$numerator_samples)]; den <- smp[split_ids(c1$denominator_samples)]
    expect_true(all(den$phase %in% c("control", "chronic")), info = c1$contrast_id)
    if (c1$contrast_class == "recovery") {
      expect_true(all(num$phase == "recovery"), info = c1$contrast_id)
    } else {
      expect_lt(max(as.numeric(num$temperature_C)), min(as.numeric(den$temperature_C)), label = c1$contrast_id)
    }
  }
})

test_that("classes follow the design: freezing below 0 C, chronic only for growth from sowing", {
  expect_setequal(unique(cons$contrast_class), c("chilling", "freezing", "recovery", "germination", "chronic"))
  expect_true(all(as.numeric(cons[contrast_class == "freezing", cold_C]) < 0))
  expect_true(all(as.numeric(cons[contrast_class == "chilling", cold_C]) >= 0))
  expect_setequal(unique(cons[contrast_class == "chronic", study_id]), "Crow2020")
})

test_that("replicated contrasts have at least two samples per side", {
  r <- cons[replicated == "TRUE"]
  expect_true(all(lengths(lapply(r$numerator_samples, split_ids)) >= 2))
  expect_true(all(lengths(lapply(r$denominator_samples, split_ids)) >= 2))
})

test_that("time-resolved contrasts compare matched times unless tagged otherwise", {
  for (i in which(nzchar(cons$time_h))) {
    c1 <- cons[i]
    if (grepl("not_time_matched", c1$tags)) next
    tn <- unique(smp[split_ids(c1$numerator_samples), time_h]); td <- unique(smp[split_ids(c1$denominator_samples), time_h])
    if (identical(td, "")) next                                 # a single untimed control, taken as matched
    expect_equal(tn, td, info = c1$contrast_id)
  }
})

test_that("display-only designs are never pooled", {
  dis <- cons[grepl("combined_treatment|secondary_control|distal|local_only|timing_unknown", tags)]
  expect_gt(nrow(dis), 0)
  expect_true(all(dis$pool == "FALSE"))
})

test_that("DESeq2 log2FC follows the numerator / denominator orientation in every kind of fit", {
  de <- pq("results", "de.parquet")
  sm <- tsv("results", "de_summary.tsv")
  ids <- c(groups = "Zhao2021_chilling_Root", groups = "Li2020a_chilling_Leaf_4C",
           blocked = "Avila2018_chilling_CG102_Leaf_24h", panel = "Zhou2022_chilling_NMpanel_Leaf_NM_25h")
  for (k in seq_along(ids)) {
    cid <- ids[[k]]
    expect_match(sm[contrast_id == cid, model], names(ids)[k], label = cid)
    d <- de[contrast_id == cid & base_mean > 50]
    r <- cor(d$log2fc, log2((d$tpm_numerator_mean + 1) / (d$tpm_denominator_mean + 1)))
    expect_gt(r, 0.8, label = cid)
  }
})
