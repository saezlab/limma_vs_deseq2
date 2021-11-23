library(edgeR)
library(readr)

GSE119931_PANC1_FOXA2KO_genes_counts <- as.data.frame(
  read_delim("data/GSE119931_PANC1.FOXA2KO.genes.counts.txt",
             "\t", escape_double = FALSE, trim_ws = TRUE))
count_df <- GSE119931_PANC1_FOXA2KO_genes_counts[,c(7:12)]
row.names(count_df) <- GSE119931_PANC1_FOXA2KO_genes_counts$Geneid

#### Pre-processing and normalisation

#First create a dataframe to summarise experimental design called targets

targets <- as.data.frame(matrix(NA,length(names(count_df)),2))
names(targets) <- c("sample","condition")
targets$sample <- names(count_df)
targets$condition <- gsub(".Rep[0-9]$","",targets$sample)
row.names(targets) <- targets$sample
targets <- targets[,2,drop = F]

y <- DGEList(counts=count_df, group=targets$condition)

keep <- filterByExpr(y, min.count = 45)

y <- y[keep, , keep.lib.sizes=FALSE]

y <- estimateDisp(y)

design <- model.matrix(~targets$condition)

fit <- glmQLFit(y, design)

qlf.2vs1 <- glmQLFTest(fit, coef=2)

top_table <- as.data.frame(topTags(qlf.2vs1, n = sum(keep)))

t <- sign(qlf.2vs1$table$logFC) * sqrt(qlf.2vs1$table$F)
# z <- zscoreT(t.stat, df=qlf.2vs1$df.total)
names(t) <- row.names(qlf.2vs1$table)
t <- data.frame(t)

top_table$gene <- row.names(top_table)
t$gene <- row.names(t)

top_table <- merge(top_table, t)

write_csv(top_table, file = "results/edgeR_DA.csv")
