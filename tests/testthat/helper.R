# Shared setup for the MaizExpressC tests (run with Rscript tests/testthat.R from the project root).
suppressPackageStartupMessages({ library(data.table); library(arrow) })
ROOT <- normalizePath(file.path(testthat::test_path(), "..", ".."))
pq <- function(...) as.data.table(read_parquet(file.path(ROOT, ...)))
tsv <- function(...) fread(file.path(ROOT, ...), colClasses = "character", na.strings = NULL)
split_ids <- function(x) setdiff(strsplit(x, ";", fixed = TRUE)[[1]], "")
GENE <- c(nii1 = "Zm00001eb255880", nnr1 = "Zm00001eb176470", gln2 = "Zm00001eb054990", nrt2 = "Zm00001eb209670")   # search tests
COLD <- c(ereb36 = "Zm00001eb318890", dbp4 = "Zm00001eb269420", ereb65 = "Zm00001eb103110", dhn3 = "Zm00001eb187010")   # cold controls (annotation/cold_positive_controls.tsv)
