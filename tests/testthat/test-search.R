# Gene ID search (HANDOFF M5): any identifier or symbol, case-insensitive, every match returned.

al <- pq("annotation", "gene_aliases.parquet")
hits <- function(q) unique(al[alias_lower == tolower(q), gene_id])

test_that("v5, v4 and v3 IDs and the symbol all find nii1", {
  for (q in c("Zm00001eb255880", "Zm00001d018161", "GRMZM2G102959", "nii1")) expect_true(GENE[["nii1"]] %in% hits(q), info = q)
})

test_that("search is case-insensitive", {
  expect_identical(hits("NII1"), hits("nii1"))
  expect_identical(hits("zm00001eb255880"), hits("Zm00001eb255880"))
})

test_that("Arabidopsis names find their maize orthologs, including dotted forms of colon names", {
  expect_true(GENE[["nii1"]] %in% hits("NIR1"))
  expect_true(GENE[["nii1"]] %in% hits("AT2G15620"))
  expect_true(GENE[["nrt2"]] %in% hits("NRT2.1"))
})

test_that("an identifier that maps to several genes returns all of them", {
  expect_gt(length(hits("NIA1")), 1)
  expect_gt(length(hits("gln1.3")), 1)
})

test_that("every quantified gene is in the annotation", {
  tables <- Sys.glob(file.path(ROOT, "results", "*", "salmon.merged.gene_counts.tsv"))
  skip_if(length(tables) == 0, "no gene tables yet")
  g <- pq("annotation", "genes.parquet")
  counts <- fread(tables[1], select = "gene_id")
  expect_true(all(counts$gene_id %in% g$gene_id))
})

test_that("cold genes are found by symbol, old IDs and Arabidopsis name", {
  for (q in c("dbp4", "Zm00001d036003", "GRMZM2G380377")) expect_true(COLD[["dbp4"]] %in% hits(q), info = q)
  expect_true(COLD[["ereb65"]] %in% hits("DDF2"))              # DREB1/CBF genes via their Arabidopsis best hit
  expect_gte(length(hits("ddf2")), 5)
})
