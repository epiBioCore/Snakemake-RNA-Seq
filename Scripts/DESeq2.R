log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(BiocParallel)
library(DESeq2)
library(tidyverse)
source(file.path(snakemake@scriptdir, 'Common.R') )

################################################################################
####################  Functions  ###############################################
################################################################################

## This Function will allow us to make a MA plot from a list of results objects



##function for adding the means of depth norm counts ,or fpkm per group the DEtable 
addMeans<-function(x,means,groups) {

  detable<-as.data.frame(x)
  
  head(means)
  groups
  ##match the gene names in the mean df with the results df
  group_means<-means[,groups]

  group_means<-group_means[match(rownames(detable),rownames(group_means)),]
  detable<-cbind(detable,group_means)
  return(detable)}


################################################################################

species_pkg <- snakemake@params[["species"]]



cat("species package")
species_pkg
snakemake@input[["species_anno"]]


load_bioconductor_package(snakemake@input[["species_anno"]], species_pkg)

##################################################################################
dds_init <- snakemake@input[["dds"]]
de<-snakemake@output[["de"]]
ma <-snakemake@output[["ma"]]
counts <- snakemake@input[["counts"]]
fpkms <- snakemake@input[["fpkms"]]
species<-snakemake@params[["org"]]
threads <- snakemake@resources[["cpus"]]
contrast <- snakemake@wildcards$contrast
alpha <- snakemake@params[["fdr"]]
ids <- snakemake@params[["ids"]]

load(dds_init)

dds
## parse the design contrast
# a contrast is always in the form 'batch+condition_group1_vs_group2', where batch(+) is optional
batch <- NA
contr <- contrast
if (grepl('\\+', contrast)) {
  batch <- strsplit(contrast, '\\+')[[1]][1]
  contr <- strsplit(contrast, '\\+')[[1]][2]
}
contr1 <- strsplit(contr, '_')[[1]]
condition <- contr1[1]
groups <- sub(paste0(condition,"_"),"",contr)
groups <- strsplit(groups,"_vs_")[[1]]

rm(contr)


colData(dds)[,"condition"] <- colData(dds)[condition]
colData(dds)$condition <- as.factor(colData(dds)$condition)
if (!is.na(batch)) {
colData(dds)[,"batch"] <- colData(dds)[batch]
colData(dds)$batch  <- as.factor(colData(dds)$batch)

}



parallel <- FALSE
if (threads > 1) {
  register(MulticoreParam(threads))
  parallel <- TRUE
}

design(dds) <-  if (!is.na(batch)){~ batch + condition} else {~ condition}
### should update to using shrinkage
dds <- DESeq(dds, parallel=parallel,betaPrior=T)

cat('batches & contrasts:\n', resultsNames(dds), '\n\n')
DE_contrast <- c("condition", groups[1], groups[2])
cat('selected contrast:', DE_contrast, '\n\n')

Res <- results(dds,contrast=DE_contrast,alpha = alpha)
summary(Res)
## fix
pdf(file = ma)
plotMA(Res, main=DE_contrast,alpha=alpha, ylim=c(-6,6))
abline(h=c(2,-2), col='red')
dev.off()



## Now let's order the result object, subset by adjusted p-value, print a 
## summary and rearrange things a bit to save these results to file
Res <- Res[order(Res$padj), ]

## add count 
count_file <- read.csv(counts,header=T,row.names=1)
Res_with_means <- addMeans(Res,count_file,groups)
colnames(Res_with_means)[7:8] <- paste(colnames(Res_with_means)[7:8],"count.mean",sep = ".")


fpkm_file <- read.csv(fpkms,header=T,row.names=1)
Res_with_means_fpkm <- addMeans(Res_with_means,fpkm_file,groups)
colnames(Res_with_means_fpkm)[9:10] <- paste(colnames(Res_with_means_fpkm)[9:10],"fpkm.mean",sep = ".")


## make the gene ids a column

Res_with_means_fpkm$GeneID <- rownames(Res_with_means_fpkm)
Res_with_means_fpkm <- dplyr::select(Res_with_means_fpkm,GeneID,everything())

Res_with_means_fpkm$GeneID <- rownames(Res_with_means_fpkm)
Res_with_means_fpkm <- dplyr::select(Res_with_means_fpkm,GeneID,everything())



if (toupper(ids) == "SYMBOL") {
  key <- "SYMBOL"
  
  columns <- c("ENSEMBL", "ENTREZID","GENENAME")  

  annot <- AnnotationDbi::select(get(species_pkg), keys=Res_with_means_fpkm$GeneID, columns=columns, keytype=key)
  annot <- filter(annot,!duplicated(SYMBOL))
  Res_to_export <- left_join(Res_with_means_fpkm,annot,by=c("GeneID" = "SYMBOL"))

} else if (toupper(ids) == "ENSEMBL") {
  key <- "ENSEMBL"
  cols <- c("SYMBOL", "ENTREZID","GENENAME")

  annot <- select(get(species_pkg), keys=Res_with_means_fpkm$GeneID, columns=cols, keytype=key)
  annot <- filter(annot,!duplicated(ENSEMBL))
  Res_to_export <- left_join(Res_with_means_fpkm,annot,by=c("GeneID" = "ENSEMBL"))

}
 




sigRes <- subset(Res_to_export, padj<=alpha)
if (nrow(sigRes)>0) {
file_name <- sub("ALL","SIG",de)
write.csv(sigRes,file=file_name,row.names = F)
 }

 ###export all results

write.csv(Res_to_export,file=de,row.names = F)







