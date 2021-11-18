library(readr)

source("scripts/support_functions.R")

DESeq2_DA <- as.data.frame(
  read_csv("results/DESeq2_DA.csv"))

vsn_limma_DA <- as.data.frame(
  read_csv("results/vsn_limma_DA.csv"))

comparison_p <- merge(DESeq2_DA[,c(1,6)], vsn_limma_DA[,c(1,5)], by = "ID")

cor.test(comparison_p$pvalue, comparison_p$P.Value, method = "spearman")

comparison_p$cohenrent_005 <- ifelse((comparison_p$pvalue < 0.05 & comparison_p$P.Value < 0.05) |
                                       (comparison_p$pvalue > 0.05 & comparison_p$P.Value > 0.05), 1, 0)
sum(comparison_p$cohenrent_005, na.rm = T)

diff_genes <- DESeq2_DA[!(DESeq2_DA$ID %in% vsn_limma_DA$ID), "ID"]
save(diff_genes,file = "results/diff_genes_limma_DESeq2.RData")

