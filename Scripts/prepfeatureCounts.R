
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

count_file<-snakemake@input[["counts"]]
counts<-read.delim(count_file,header=TRUE,stringsAsFactors=F,skip=1,check.names=F)

head(counts)
rownames(counts)<-counts$Geneid
annotation <- counts[,2:6]
counts<-counts[,-c(1:6)]

colnames(counts)<-gsub("_sorted.bam","",basename(colnames(counts)))
head(counts)

write.csv(counts, file = snakemake@output[["counts"]])
write.csv(annotation,file = snakemake@output[["annot"]])
