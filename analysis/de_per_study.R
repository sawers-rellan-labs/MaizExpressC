#!/usr/bin/env Rscript
# M4: per-study differential expression for every contrast in metadata/contrasts.tsv.
#
#   Rscript analysis/de_per_study.R              # all studies with gene tables in results/
#   Rscript analysis/de_per_study.R Zhao2021 ... # selected studies (writes results/de_<study>.parquet only)
#
# As in MaizExpressN: DESeq2 per study, design ~ 0 + grp, where grp is the sample group (the
# same group key build_contrasts.py uses) and the model holds every included sample of the
# study. Libraries in metadata/library_exclusions.tsv are left out. Each contrast is a Wald test
# of the mean of its numerator groups against the mean of its denominator groups, so
# log2FC = log2(cold / control) (recovery: recovered / control), as in contrasts.tsv. Counts
# come from salmon via tximport with average-transcript-length offsets; 3' tag studies
# (Xue2021) use counts without length offsets.
# Two cold designs get their own fit (column `model` in de_summary.tsv):
#   blocked: contrasts with a batch_term (Avila2018, Makarevitch2015, Zhang2017's time course)
#     are fitted on their own libraries with ~ 0 + grp + batch, batch = replicate (separate
#     plantings and sequencing runs); without the block if it is not estimable after exclusions.
#   panel: Zhou2022's diversity-panel contrasts (tag panel) compare cold with control across
#     inbreds with one library each, paired by inbred: ~ genotype + cond on those libraries
#     (inbreds whose partner library was excluded are left out).
# Unshrunk log2FC and lfcSE are kept for pooling; ashr-shrunk log2FC is for display (apeglm
# cannot shrink arbitrary contrasts). Unreplicated contrasts get the log2 ratio of mean
# normalised counts (+0.5) and no SE, p or padj.
# single_library (replicated contrasts; MaizExpressN rule, Rubén 22 Oct 2026): the fold change
# rests on one library. TRUE when one library holds more than 90% of the gene's normalised counts
# across the contrast's libraries, or, for a gene with padj < 0.05, when the higher condition
# without its top library is no higher than the lower condition (mean normalised counts). Like
# Cook's outliers, flagged genes get pvalue and padj NA (other genes' padj is unchanged); log2FC,
# SE and shrunk log2FC stay for display. Pooling leaves out every replicated row without a p-value.
#
# Outputs: results/de.parquet (gene x contrast) and results/de_summary.tsv (one row per
# contrast: model, samples used, genes tested, DEG counts, single_library genes).

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(arrow)
  library(BiocParallel)
})
stopifnot(requireNamespace("ashr", quietly = TRUE))  # used by lfcShrink(type = "ashr")
options(warn = 1)
register(MulticoreParam(as.integer(Sys.getenv("DE_WORKERS", "6"))))   # DESeq(parallel = TRUE) splits genes

NO_LENGTH <- c("Xue2021")         # 3' tag-seq (QuantSeq): no length offsets
PADJ_DEG <- 0.05; LFC_DEG <- 1      # DEG call for the summary table only
MIN_COUNT <- 10                    # a gene is tested if >= MIN_COUNT reads in >= MIN_SAMPLES samples
SPIKE_SHARE <- 0.9; SPIKE_PADJ <- 0.05   # single_library rule (see header)
GROUP_KEY <- c("study_id", "genotype", "tissue_group", "tissue_detail", "temp_label_paper", "temperature_C",
               "temperature_night_C", "control_temperature_C", "phase", "time_h", "recovery_h", "light_condition",
               "other_factor", "other_factor_level")

samples <- fread("metadata/samples.tsv", colClasses = "character", na.strings = NULL)
cons <- fread("metadata/contrasts.tsv", colClasses = "character", na.strings = NULL)
excl <- fread("metadata/library_exclusions.tsv", colClasses = "character")
samples[, group := do.call(paste, c(.SD, sep = "|")), .SDcols = GROUP_KEY]
use <- samples[include == "TRUE" & !(sample_id %in% excl$sample_id)]

split_ids <- function(x) setdiff(strsplit(x, ";", fixed = TRUE)[[1]], "")
has_tag <- function(tags, tag) grepl(paste0("(^|;)", tag, "(;|$)"), tags)

# N, D: normalised counts (genes x libraries) of the numerator and denominator; see the header
is_single_library <- function(N, D, lfc, padj) {
  sn <- rowSums(N); sd_ <- rowSums(D); mn <- apply(N, 1, max); md <- apply(D, 1, max)
  share <- pmax(mn, md) / (sn + sd_)
  up <- lfc > 0
  rest_hi <- ifelse(up, (sn - mn) / (ncol(N) - 1), (sd_ - md) / (ncol(D) - 1))   # higher condition without its top library
  lo <- ifelse(up, sd_ / ncol(D), sn / ncol(N))
  flag <- (sn + sd_ > 0 & share > SPIKE_SHARE) | (!is.na(padj) & padj < SPIKE_PADJ & rest_hi <= lo)
  flag %in% TRUE
}

# gene x sample matrices from each owning study's nf-core tables, cached
FILES <- c(counts = "salmon.merged.gene_counts.tsv", lengths = "salmon.merged.gene_lengths.tsv",
           tpm = "salmon.merged.gene_tpm.tsv")
cache <- new.env()
read_mat <- function(owner, kind) {
  key <- paste(owner, kind)
  if (is.null(cache[[key]])) {
    d <- fread(file.path("results", owner, FILES[[kind]]))
    m <- as.matrix(d[, -(1:2)]); rownames(m) <- d$gene_id
    cache[[key]] <- m
  }
  cache[[key]]
}
get_mat <- function(sids, kind) {
  owners <- samples[match(sids, sample_id), study_id]
  parts <- lapply(unique(owners), function(o) read_mat(o, kind)[, sids[owners == o], drop = FALSE])
  genes <- rownames(parts[[1]])
  m <- do.call(cbind, lapply(parts, function(p) p[genes, , drop = FALSE]))
  m[, sids, drop = FALSE]
}

build_dds <- function(st, meta, design) {
  sids <- rownames(meta)
  cnt <- get_mat(sids, "counts")
  if (st %in% NO_LENGTH) return(DESeqDataSetFromMatrix(round(cnt), meta, design))
  len <- get_mat(sids, "lengths")
  len[!is.finite(len) | len <= 0] <- NA
  rm <- rowMeans(len, na.rm = TRUE); rm[!is.finite(rm)] <- 1
  idx <- which(is.na(len), arr.ind = TRUE); len[idx] <- rm[idx[, 1]]
  txi <- list(counts = cnt, abundance = get_mat(sids, "tpm"), length = len, countsFromAbundance = "no")
  suppressMessages(DESeqDataSetFromTximport(txi, meta, design))
}

# one DESeq2 fit: kind "groups" (~ 0 + grp), "blocked" (~ 0 + grp + batch) or "panel" (~ genotype + cond).
# Groups with a single library are left out of the GLM fit: in ~ 0 + grp a one-library group fits its
# library exactly and adds no residual df, but each one adds a coefficient (Zhou2022 has 143 of 173).
# Gene filter and normalised counts for unreplicated contrasts use all the fit's libraries, as before.
fit_model <- function(st, sids, kind, cold = NULL) {
  meta <- as.data.frame(use[match(sids, sample_id)])
  rownames(meta) <- sids
  note <- ""
  if (kind == "panel") {
    meta$cond <- factor(ifelse(sids %in% cold, "cold", "control"), levels = c("control", "cold"))
    meta$genotype <- factor(meta$genotype)
    design <- ~ genotype + cond
    infit <- rep(TRUE, nrow(meta))
  } else {
    lev <- unique(meta$group)
    meta$grp <- factor(paste0("L", match(meta$group, lev)), levels = paste0("L", seq_along(lev)))
    meta$batch <- factor(meta$replicate)
    design <- ~ 0 + grp
    infit <- meta$grp %in% names(which(table(meta$grp) >= 2))
  }
  fmeta <- droplevels(meta[infit, , drop = FALSE])
  if (kind != "panel" && nlevels(fmeta$grp) < 2) fmeta <- fmeta[0, , drop = FALSE]   # one group: nothing to test
  if (kind == "blocked" && nrow(fmeta)) {
    mm <- model.matrix(~ 0 + grp + batch, fmeta)
    if (nlevels(fmeta$batch) > 1 && qr(mm)$rank == ncol(mm)) design <- ~ 0 + grp + batch
    else note <- "replicate block not estimable after exclusions: fitted without it"
  }
  replicated_fit <- nrow(fmeta) > 0 && nrow(fmeta) > qr(model.matrix(design, fmeta))$rank
  all <- build_dds(st, meta, ~ 1)
  keep <- rowSums(counts(all) >= MIN_COUNT) >= (if (replicated_fit) 2 else 1)
  all <- estimateSizeFactors(all[keep, ])
  dds <- if (replicated_fit) DESeq(build_dds(st, fmeta, design)[keep, ], quiet = TRUE, parallel = TRUE)
  list(dds = dds, ncounts = if (replicated_fit) counts(dds, normalized = TRUE), ncounts_all = counts(all, normalized = TRUE),
       tpm = get_mat(sids, "tpm"), kind = kind,
       gmap = if (kind == "panel") NULL else setNames(as.character(fmeta$grp), rownames(fmeta)),
       replicated_fit = replicated_fit, note = note)
}

run_study <- function(st) {
  cs <- cons[study_id == st]
  missing <- st[!file.exists(file.path("results", st, FILES[["counts"]]))]
  if (length(missing)) return(list(note = paste("no gene table for", missing)))
  cs[, fit_key := fifelse(has_tag(tags, "panel"), paste0("panel:", contrast_id),
                          fifelse(batch_term != "", paste0("blocked:", batch_term), "groups"))]
  study_sids <- use[study_id == st, sample_id]
  fits <- list(); out <- list(); summ <- list()
  for (i in seq_len(nrow(cs))) {
    c1 <- cs[i]
    num <- intersect(split_ids(c1$numerator_samples), use$sample_id)
    den <- intersect(split_ids(c1$denominator_samples), use$sample_id)
    dropped <- setdiff(c(split_ids(c1$numerator_samples), split_ids(c1$denominator_samples)), c(num, den))
    if (startsWith(c1$fit_key, "panel:")) {          # keep inbreds with both libraries
      g <- setNames(use$genotype, use$sample_id)
      both <- intersect(g[num], g[den])
      dropped <- c(dropped, num[!g[num] %in% both], den[!g[den] %in% both])
      num <- num[g[num] %in% both]; den <- den[g[den] %in% both]
    }
    if (!length(num) || !length(den)) {
      summ[[i]] <- data.table(contrast_id = c1$contrast_id, study_id = st, note = "no samples left after exclusions")
      next
    }
    key <- c1$fit_key
    if (is.null(fits[[key]])) {
      kind <- sub(":.*", "", key)
      sids <- switch(kind,
                     groups = unique(c(study_sids, unlist(lapply(c(cs[fit_key == key, numerator_samples],
                                                                   cs[fit_key == key, denominator_samples]), split_ids)))),
                     blocked = unlist(lapply(c(cs[fit_key == key, numerator_samples],
                                               cs[fit_key == key, denominator_samples]), split_ids)),
                     panel = c(num, den))
      sids <- unique(sids[sids %in% use$sample_id])
      fits[[key]] <- fit_model(st, sids, kind, cold = num)
    }
    f <- fits[[key]]; dds <- f$dds; ncounts <- f$ncounts
    rep_ok <- c1$replicated == "TRUE" && f$replicated_fit && length(num) >= 2 && length(den) >= 2
    if (rep_ok) {
      if (f$kind == "panel") {
        res <- results(dds, name = "cond_cold_vs_control", alpha = PADJ_DEG)
      } else {
        ng <- unique(f$gmap[num]); dg <- unique(f$gmap[den])
        vec <- setNames(numeric(length(resultsNames(dds))), resultsNames(dds))
        vec[paste0("grp", ng)] <- 1 / length(ng); vec[paste0("grp", dg)] <- -1 / length(dg)
        res <- results(dds, contrast = vec, alpha = PADJ_DEG)
      }
      sh <- suppressMessages(lfcShrink(dds, res = res, type = "ashr", quiet = TRUE))
      d <- data.table(gene_id = rownames(res), base_mean = res$baseMean, log2fc = res$log2FoldChange,
                      lfc_se = res$lfcSE, log2fc_shrunk = sh$log2FoldChange, stat = res$stat,
                      pvalue = res$pvalue, padj = res$padj)
      d[, single_library := is_single_library(ncounts[gene_id, num, drop = FALSE], ncounts[gene_id, den, drop = FALSE],
                                           log2fc, padj)]
      d[single_library == TRUE, `:=`(pvalue = NA_real_, padj = NA_real_)]
    } else {
      ncounts <- f$ncounts_all
      a <- rowMeans(ncounts[, num, drop = FALSE]); b <- rowMeans(ncounts[, den, drop = FALSE])
      d <- data.table(gene_id = rownames(ncounts), base_mean = rowMeans(ncounts[, c(num, den), drop = FALSE]),
                      log2fc = log2((a + 0.5) / (b + 0.5)), lfc_se = NA_real_, log2fc_shrunk = NA_real_,
                      stat = NA_real_, pvalue = NA_real_, padj = NA_real_, single_library = NA)
    }
    tp <- f$tpm[d$gene_id, , drop = FALSE]
    d[, `:=`(contrast_id = c1$contrast_id, study_id = st,
             tpm_numerator_mean = rowMeans(tp[, num, drop = FALSE]),
             tpm_denominator_mean = rowMeans(tp[, den, drop = FALSE]),
             replicated = rep_ok)]
    out[[i]] <- d
    summ[[i]] <- data.table(
      contrast_id = c1$contrast_id, study_id = st, contrast_class = c1$contrast_class,
      model = paste(f$kind, if (f$kind == "blocked" && !nzchar(f$note)) "(replicate)" else ""),
      replicated = rep_ok, n_numerator = length(num), n_denominator = length(den),
      excluded_libraries = paste(dropped, collapse = ";"), genes_tested = nrow(d),
      n_up = if (rep_ok) d[padj < PADJ_DEG & log2fc >= LFC_DEG, .N] else NA_integer_,
      n_down = if (rep_ok) d[padj < PADJ_DEG & log2fc <= -LFC_DEG, .N] else NA_integer_,
      n_single_library = if (rep_ok) d[single_library == TRUE, .N] else NA_integer_,
      note = trimws(paste(f$note, if (f$replicated_fit) "" else "no replicates in this fit: log2 ratio of normalised counts")))
  }
  list(de = rbindlist(out), summary = rbindlist(summ, fill = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
studies <- if (length(args)) args else sort(unique(cons$study_id))
res <- list(); summ <- list()
for (st in studies) {
  t0 <- Sys.time()
  r <- tryCatch(run_study(st), error = function(e) list(note = paste("error:", conditionMessage(e))))
  if (!is.null(r$note)) {
    message(sprintf("%-20s skipped: %s", st, r$note))
    summ[[st]] <- data.table(study_id = st, note = r$note)
    next
  }
  res[[st]] <- r$de; summ[[st]] <- r$summary
  message(sprintf("%-20s %3d contrasts  %6.1f s", st, nrow(r$summary), as.numeric(Sys.time() - t0, units = "secs")))
}
de <- rbindlist(res, use.names = TRUE)
summary <- rbindlist(summ, fill = TRUE)
setcolorder(de, c("gene_id", "contrast_id", "study_id", "base_mean", "log2fc", "lfc_se", "log2fc_shrunk",
                  "stat", "pvalue", "padj", "single_library", "tpm_numerator_mean", "tpm_denominator_mean", "replicated"))
if (length(args)) {
  for (st in names(res)) {
    write_parquet(res[[st]], file.path("results", paste0("de_", st, ".parquet")))
    fwrite(summ[[st]], file.path("results", paste0("de_summary_", st, ".tsv")), sep = "\t")
  }
} else {
  write_parquet(de, "results/de.parquet")
  fwrite(summary, "results/de_summary.tsv", sep = "\t")
}
message(sprintf("%d gene x contrast rows from %d studies", nrow(de), length(res)))
