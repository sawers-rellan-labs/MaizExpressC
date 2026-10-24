#!/usr/bin/env Rscript
# M5: cold positive controls -> annotation/cold_positive_controls.tsv
#
#   Rscript annotation/build_cold_controls.R      (after analysis/de_per_study.R)
#
# The gene list below is curated from the annotation, not from the data: maize DREB1/CBF genes
# (MaizeGDB AP2-EREBP 'ereb' / 'dbp' / 'dbf' genes whose Arabidopsis best hit is DDF2/DREB1F, a CBF
# subfamily member; Ensembl lists no 1:1 maize orthologue of CBF1-4; ereb23 has no Ensembl hit, but its
# EnTAP best match is Setaria DREB1A and it lies 7 kb from ereb36 in a tandem cluster), dehydrins (COR/LEA genes),
# and the three SCREAM-like bHLH genes as ICE1 homologues. Each gene's role says what is expected:
#   positive_control: up in leaf and shoot chilling or freezing at <= 10 C (the DREB1/CBF genes);
#   weak_control: up in many of those contrasts, smaller fold changes (dehydrins);
#   reference: ICE1 homologues, regulated after translation, so little change is expected.
# The M4 columns then record how each gene behaved in the replicated, pooled contrasts of
# results/de.parquet: aerial_acute = class chilling or freezing, tissue Leaf, Shoot or Whole
# seedling, cold_C <= 10; root = chilling in roots; chronic = Crow2020. up / down = padj < 0.05.
# DREB1 genes are barely expressed in roots and do not respond to chronic 22 vs 32 C growth.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

CONTROLS <- rbindlist(list(
  data.table(group = "DREB1/CBF", role = "positive_control",
             gene_id = c("Zm00001eb103110", "Zm00001eb318890", "Zm00001eb269420", "Zm00001eb318900",
                         "Zm00001eb250560", "Zm00001eb073550", "Zm00001eb037530", "Zm00001eb103100")),
  data.table(group = "dehydrin (COR/LEA)", role = "weak_control",
             gene_id = c("Zm00001eb187010", "Zm00001eb376710", "Zm00001eb052790", "Zm00001eb250120")),
  data.table(group = "ICE1 homologue (SCREAM-like bHLH)", role = "reference",
             gene_id = c("Zm00001eb113220", "Zm00001eb142500", "Zm00001eb170530"))
))

genes <- as.data.table(read_parquet("annotation/genes.parquet"))
cons <- fread("metadata/contrasts.tsv", colClasses = "character", na.strings = NULL)
de <- as.data.table(read_parquet("results/de.parquet", col_select = c("gene_id", "contrast_id", "log2fc", "padj", "replicated")))
de <- de[gene_id %in% CONTROLS$gene_id & replicated == TRUE]
de <- merge(de, cons[pool == "TRUE", .(contrast_id, contrast_class, tissue_group, cold_C = as.numeric(cold_C))], by = "contrast_id")
de[, set := fifelse(contrast_class %in% c("chilling", "freezing") & tissue_group %in% c("Leaf", "Shoot", "Whole seedling") &
                      cold_C <= 10, "aerial_acute",
           fifelse(contrast_class == "chilling" & tissue_group == "Root", "root",
           fifelse(contrast_class == "chronic", "chronic", NA_character_)))]
s <- de[!is.na(set), .(tested = sum(!is.na(log2fc)), up = sum(padj < 0.05 & log2fc > 0, na.rm = TRUE),
                       down = sum(padj < 0.05 & log2fc < 0, na.rm = TRUE), median_log2fc = round(median(log2fc, na.rm = TRUE), 2)),
        by = .(gene_id, set)]
w <- dcast(s, gene_id ~ set, value.var = c("tested", "up", "down", "median_log2fc"))
out <- merge(CONTROLS, genes[, .(gene_id, symbol, v4_ids, description, ath_best_hit, ath_symbol)], by = "gene_id", sort = FALSE)
out <- merge(out, w, by = "gene_id", all.x = TRUE, sort = FALSE)
cols <- c("gene_id", "symbol", "group", "role", "v4_ids", "description", "ath_best_hit", "ath_symbol",
          as.vector(outer(c("tested_", "up_", "down_", "median_log2fc_"), c("aerial_acute", "root", "chronic"), paste0)))
setcolorder(out, intersect(cols, names(out)))
out <- out[match(CONTROLS$gene_id, gene_id)]
fwrite(out, "annotation/cold_positive_controls.tsv", sep = "\t")
n_acute <- uniqueN(de[set == "aerial_acute", contrast_id])
message(sprintf("%d genes -> annotation/cold_positive_controls.tsv (%d aerial acute contrasts)", nrow(out), n_acute))
print(out[, .(symbol, role, acute = sprintf("%d up / %d down of %d, median %+.1f", up_aerial_acute, down_aerial_acute,
                                             tested_aerial_acute, median_log2fc_aerial_acute),
              root = sprintf("%d / %d of %d", up_root, down_root, tested_root))])
