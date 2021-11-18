library(dorothea)
library(decoupleR)
library(readr)
library(org.Hs.eg.db)

min_reg_size <- 1

DESeq2_DA <- as.data.frame(
  read_csv("results/DESeq2_DA.csv"))
names(DESeq2_DA)[5] <- "t"
DESeq2_DA$t <- DESeq2_DA$t * -1 #Pau's 'ol fliperoo

vsn_limma_DA <- as.data.frame(
  read_csv("results/vsn_limma_DA.csv"))

vsn_limma_DA_with_low_expr <- as.data.frame(
  read_csv("results/vsn_limma_DA_with_low_expr.csv"))

ttops <- list(DESeq2_DA, vsn_limma_DA, vsn_limma_DA_with_low_expr)

dorothea <- dorothea_hs
dorothea <- dorothea[dorothea$confidence %in% c("A","B","C"),]
dorothea$likelihood <- 1

TF_res <- list()
inputs <- list()
i <- 1
for(ttop in ttops)
{
  ttop$ID <- gsub("[.].*","",ttop$ID)
  mapping <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ttop$ID, 'SYMBOL', 'ENSEMBL')
  ttop$ID <- mapping[ttop$ID]

  dubs <- ttop[duplicated(ttop$ID),"ID"]

  ttop <- ttop[!ttop$ID %in% dubs,]

  row.names(ttop) <- ttop $ID
  ttop <- ttop[,"t",drop = F]

  inputs[[i]] <- ttop

  regulons <- intersect_regulons(mat = ttop, network = dorothea, .source = tf, .target = target, minsize = min_reg_size)

  TF_res[[i]] <- run_mlm(mat = ttop, network = regulons, .source = "tf")
  i <- i+1
}

TF_merge_mlm <- merge(TF_res[[1]][,c(2,4)],TF_res[[2]][,c(2,4)], by = "source")
TF_merge_mlm <- merge(TF_merge_mlm,TF_res[[3]][,c(2,4)], by = "source")
names(TF_merge_mlm) <- c("tf","DESeq2_mlm","limma_filter_mlm","limma_full_mlm")

TF_res <- list()
inputs <- list()
i <- 1
for(ttop in ttops)
{
  ttop$ID <- gsub("[.].*","",ttop$ID)
  mapping <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ttop$ID, 'SYMBOL', 'ENSEMBL')
  ttop$ID <- mapping[ttop$ID]

  dubs <- ttop[duplicated(ttop$ID),"ID"]

  ttop <- ttop[!ttop$ID %in% dubs,]

  row.names(ttop) <- ttop $ID
  ttop <- ttop[,"t",drop = F]

  inputs[[i]] <- ttop
  regulons <- intersect_regulons(mat = ttop, network = dorothea, .source = tf, .target = target, minsize = min_reg_size)

  TF_res[[i]] <- run_wmean(mat = ttop, network = regulons, .source = "tf", times = 1000)
  TF_res[[i]] <- TF_res[[i]][TF_res[[i]]$statistic == "norm_wmean",]
  i <- i+1
}

TF_merge_wmean <- merge(TF_res[[1]][,c(2,4)],TF_res[[2]][,c(2,4)], by = "source")
TF_merge_wmean <- merge(TF_merge_wmean,TF_res[[3]][,c(2,4)], by = "source")
names(TF_merge_wmean) <- c("tf","DESeq2_wmean","limma_filter_wmean","limma_full_wmean")

full_merge <- merge(TF_merge_mlm, TF_merge_wmean)
row.names(full_merge) <- full_merge$tf
full_merge <- full_merge[,-1]

pheatmap::pheatmap(full_merge, display_numbers = F, show_rownames = F)
