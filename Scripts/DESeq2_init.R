log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(BiocParallel)
library(tidyverse)


counts <- snakemake@input[["counts"]]
len <- snakemake@input[["len"]]
meta <- snakemake@params[["meta"]]
out <- snakemake@output[["dds"]]
comparisons <- snakemake@params$comparisons
replicates <- snakemake@params$replicates
threads <- snakemake@resources$cpus

if (snakemake@config$use_cuffnorm_fpkms) {fpkms <- snakemake@input[["cuff_fpkms"]]}
counts<-read.csv(counts,header=TRUE,stringsAsFactors=F,check.names = F,row.names=1)

head(counts)
samples <- read.delim(meta, sep = "\t", na.strings = "", comment.char = "#", stringsAsFactors = F)
if ("replicate" %in% colnames(samples) & isTRUE(replicates)) {
  samples$replicate[is.na(samples$replicate)] <- as.character(samples$sample[is.na(samples$replicate)])
  samples <- subset(samples, !duplicated(replicate))
  row.names(samples) <- samples$replicate
} else {
  row.names(samples) <- samples$sample
}
## remove any white spaces
samples[,1:ncol(samples)] <- lapply(samples[,1:ncol(samples)],trimws) 
colData <- samples

parallel <- FALSE
if (threads > 1) {
  register(MulticoreParam(threads))
  parallel <- TRUE
}

##match rownames of colData with the colnames of count matrix
##check if rownames and colnames match

counts <- counts[,match(rownames(colData),colnames(counts))]
## Create DESeqDataSet
## Consolidates the data into an object that can be utilized by the program for 
## calculations and plotting.  Requires you to specify the experimental design 
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=colData,
                              design= ~ 1)
## ~ 1 this is initial dds only, will reset design when using for each contrast
## add gene lengths to dds

gene_lengths <- read.csv(len,header=T,row.names=1)
head(gene_lengths)
gene_lengths <- gene_lengths[match(rownames(counts),rownames(gene_lengths)),]
colnames(gene_lengths)[5] <- "basepairs"
head(gene_lengths)
mcols(dds) <- DataFrame(gene_lengths)
dds <- DESeq(dds,parallel=parallel)
save(dds,file=out)

####################  Count and FPKM data   ####################################

### get a list of groups 

c_sheet <- read.delim(comparisons,header=T,na.strings=c(""," "))
## comparisons file should have 4 columns:
## Batch/covariate Column Group1 Group2
## where Column is the column name in sample sheet containing groups 1 and 2 
## or 3 columns Column Group1 Group2
head(c_sheet)
contrasts <- character(0)
if (ncol(c_sheet) == 4) {
  for (i in 1:nrow(c_sheet)) {
      batch <- trimws(c_sheet[i,1],which = 'both')
      col <- trimws(c_sheet[i,2],which = 'both')
      g1  <- trimws(c_sheet[i,3],which = 'both')
      g2  <- trimws(c_sheet[i,4],which = 'both')
      ## if contrasts with batch and without batch are in same sheet
      if (is.na(batch)) {
          contrast <- sprintf("%s_%s_vs_%s",col,g1,g2)
          cat(contrast)
      } else {
          contrast <- sprintf("%s+%s_%s_vs_%s", batch,col,g1,g2)
          cat(contrast)

      } 

      contrasts <- c(contrasts,contrast)        
  }
} else if (ncol(c_sheet) == 3) {
  for (i in 1:nrow(c_sheet)) {
    col <- trimws(c_sheet[i,1],which = 'both')
    g1  <- trimws(c_sheet[i,2],which = 'both')
    g2  <- trimws(c_sheet[i,3],which = 'both')

    contrast <- sprintf("%s_%s_vs_%s",col,g1,g2)
    contrasts <- c(contrasts,contrast)        

  }
}


contrasts
## Retrieve count data and clean up the data frame 
count_data <- counts(dds, normalized=T)
head(count_data)

if ("descriptive_name" %in% colnames(colData)) colnames(count_data) <- colData$descriptive_name
## get mean counts for each var

if (ncol(c_sheet) == 4) {
  column <- c_sheet[,2]
} else if (ncol(c_sheet) == 3) {
    column <- c_sheet[,1]
}

column 
cols <- unique(column)
cols
colData(dds)[,cols]
if (length(cols) ==1) {
count_means<-t(apply(count_data,1,function(x) tapply(x,colData(dds)[,cols],mean,na.rm=T)))
} else if (length(cols)>1) {
count_means<-lapply(cols,function(c) t(apply(count_data,1,function(x) tapply(x,colData(dds)[,cols],mean,na.rm=T))))
count_means <- do.call(cbind,count_means)
}



##count means and counts should be in the same order, but just in case....
count_means<-count_means[match(rownames(count_data),rownames(count_means)),]

count_data_with_means <- data.frame(cbind(count_data,count_means),check.names=T)
head(count_data_with_means)
write.csv(count_data_with_means, file= snakemake@output[["counts"]], quote=F)


if (!snakemake@config$use_cuffnorm_fpkms) {

fpkm_data <- fpkm(dds, robust=T)
if ("descriptive_name" %in% colnames(colData)) colnames(fpkm_data) <- colData$descriptive_name


if (length(cols) ==1) {
fpkm_means<-t(apply(fpkm_data,1,function(x) tapply(x,colData(dds)[,cols],mean,na.rm=T)))
} else if (length(cols)>1) {
fpkm_means<-lapply(cols,function(c) t(apply(fpkm_data,1,function(x) tapply(x,colData(dds)[,cols],mean,na.rm=T))))
fpkm_means <- do.call(cbind,fpkm_means)
}


##fpkm means and counts should be in the same order, but just in case....
fpkm_means<-fpkm_means[match(rownames(fpkm_data),rownames(fpkm_means)),]

fpkm_data_with_means <- data.frame(cbind(fpkm_data,fpkm_means))
head(fpkm_data_with_means)
write.csv(fpkm_data_with_means, file=snakemake@output[["fpkms"]], quote=F)

} else {
fpkms
fpkm_values <- read.table(fpkms,header=T,stringsAsFactors=F)
idx <- fpkm_values$tracking_id
fpkm_values <- as.matrix(fpkm_values[,-1])
rownames(fpkm_values) <- idx

##sample table from cuffnorm
samples <- read.delim(file.path(dirname(fpkms),"samples.table"),stringsAsFactors=F)
colnames(fpkm_values) <- gsub("_sorted.bam","",basename(samples[match(colnames(fpkm_values),samples$sample_id),"file"]))

##first match the gene names in dds object to cufflinks fpkm table
fpkm_data<-fpkm_values[match(rownames(count_data),rownames(fpkm_values)),]


fpkm_data[is.na(fpkm_data)] <- 0

##check that fpkm table is in same order as colData
fpkm_data <- fpkm_data[,match(rownames(colData(dds)),colnames(fpkm_data))]

if ("descriptive_name" %in% colnames(colData)) colnames(fpkm_data) <- colData$descriptive_name


if (length(cols) ==1) {
fpkm_means<-t(apply(fpkm_data,1,function(x) tapply(x,colData(dds)[,cols],mean,na.rm=T)))
} else if (length(cols)>1) {
fpkm_means<-lapply(cols,function(c) t(apply(fpkm_data,1,function(x) tapply(x,colData(dds)[,cols],mean,na.rm=T))))
fpkm_means <- do.call(cbind,fpkm_means)
}



#fpkm means and counts should be in the same order, but just in case....
fpkm_means<-fpkm_means[match(rownames(fpkm_data),rownames(fpkm_means)),]

fpkm_data_with_means <- data.frame(cbind(fpkm_data,fpkm_means),check.names =F)

write.csv(fpkm_data_with_means, file= snakemake@output[["fpkms"]], quote=F)

}