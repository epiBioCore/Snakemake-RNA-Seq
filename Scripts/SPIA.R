log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(DOSE)

library(SPIA)
source(file.path(snakemake@scriptdir, 'Common.R') )
species <- snakemake@params[["species"]]

#### get Spia organism
org <- substr(species,1,3)



sig_genes <- snakemake@input[["sig"]]
out <- snakemake@output[["object"]]
dir <- snakemake@params[["out"]]
comp <- snakemake@wildcards$contrast
species_pkg <- snakemake@params[["species_anno"]]
fdr <- snakemake@params[["fdr"]]

load_bioconductor_package(snakemake@input[["species_anno"]], species_pkg)

ALL_genes <- read.csv(gsub("SIG","ALL",sig_genes))
sig_genes <- read.csv(sig_genes)
##
top <- ALL_genes[!is.na(ALL_genes$ENTREZID),]
top<- top[!duplicated(top$ENTREZID),]

tgl <- subset(top,padj < fdr)
## what to do if no DE genes?
DE <-tgl[,"log2FoldChange"]
names(DE)  <- tgl$ENTREZID

ALL <- top$ENTREZID



res <- spia(de=DE,all=ALL,organism=org,plots = F,beta = NULL,combine = "fisher",verbose = "FALSE")



write.table(res,file = out,sep = "\t",row.names = F,quote = F)

