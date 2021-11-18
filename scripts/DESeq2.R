library(DESeq2)
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


dds <- DESeqDataSetFromMatrix(countData = count_df,
                              colData = targets,
                              design = ~ condition)
dds

keep <- rowSums(counts(dds)) >= 9
dds <- dds[keep,]

# dds$condition <- factor(dds$condition, levels = c("WT","KO"))

dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
res$ID <- row.names(res)
res <- res[,c(7,1,2,3,4,5,6)]

write_csv(res, file = "results/DESeq2_DA.csv")
