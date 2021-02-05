log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)

de_files <- unlist(snakemake@input)

fdr <- snakemake@params[["fdr"]]
de <- lapply(de_files,read.csv,header = T)
names(de) <- gsub("_DEtable_ALL_genes.csv","",basename(de_files))

de_summary <- data.frame(
    Downregulated = sapply(de,function(x) nrow(subset(x,padj < fdr & log2FoldChange <0))),
    Upregulated = sapply(de,function(x)nrow(subset(x,padj < fdr & log2FoldChange > 0))),
    Total = sapply(de,function(x) nrow(subset(x,padj < fdr)))
)

de_summary$Comparison<-names(de)
de_summary <- select(de_summary,Comparison,everything())

write.csv(de_summary,file=snakemake@output[["s"]],row.names=F)
