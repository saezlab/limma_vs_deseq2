#Main libraries
library(readr)
library(vsn)
library(limma)

source("scripts/support_functions.R")

### Import the raw count dataframe
#downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119931
#download the file : GSE119931_PANC1.FOXA2KO.genes.counts.txt.gz and decompress it in the data folder

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

#Make some plots to check what the data looks like after only a log2 transformation

#First we remove rows that contain only 0
count_df <- count_df[rowSums(count_df) > 0,]
#remaining 0 have to be made as NA so that log2 transformation is possible
count_df[count_df == 0] <- NA

plot(density(as.numeric(as.matrix(log2(count_df))), na.rm = T))
abline(v = 5.5)
#from the density plot we define the minimum cutoff of expression.
#Distributions are weird under log2(count) of around 4

count_df[log2(count_df) < 5.5 ] <- NA

#remove rows that don't have enough well measured genes in enough samples

count_df <- count_df[rowSums(is.na(count_df[,c(1:3)])) < 2,]
count_df <- count_df[rowSums(is.na(count_df[,c(4:6)])) < 2,]

### VSN normalisation
#now we can normalise the cleaned dataframe using vsn
fit <- vsnMatrix(as.matrix(count_df)) #train vsn parameters

#make sure the mean/sd trend is not going crazy
meanSdPlot(fit)

#if good, normalise data with the trained parameters of vsn
count_df_vsn <- as.data.frame(vsn::predict(fit,as.matrix(count_df)))

### Identifier kung-fu (optional)
#since here with have ensembl id but most our ressources are based on either uniprot or gene symbole
#we need to do some identifer kung-fu

### LIMMA differential analysis
#now let's run a simple differential analysis using a simple wrapper for such situation

#first check the conditions order
unique(targets$condition)

#we want to compare the KO with the WT so we build a comparison list
comparisons <- list(c(2,-1)) #each vector of the list represent the contrasts, here we substract the first condition (-1) to the second one (2)

#now that the comparisons are defined, we can run limma
limmaRes <- runLimma(measurements = count_df_vsn,
                     targets = targets,
                     comparisons = comparisons)

#once limma has run, we extract the statistic dataframe summarise the differential analysis
ttop_KOvsWT <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(count_df_vsn[,1]), adjust.method = "fdr"))

volcano_nice(ttop_KOvsWT, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 10, label = T)


write_csv(ttop_KOvsWT, file = "results/vsn_limma_DA.csv")
