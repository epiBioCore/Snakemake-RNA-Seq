log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(gplots)
library(RColorBrewer)

load(snakemake@input[["dds"]])
rowsl <- snakemake@params[["row_labels"]]
colsl <- snakemake@params[["col_labels"]]
rowsl
colsl
## Perform rLog transformation on the data 
rld <- rlog(dds)




hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)

##sample hm
dists <- dist(t(assay(rld)))
mat <- as.matrix(dists)
rownames(mat) <- as.character(colData(dds)[,rowsl])
colnames(mat) <- as.character(colData(dds)[,colsl])

hc <- hclust(dists)

## Plot HM and save to file.
pdf(file=snakemake@output[["hm"]])
hm <- heatmap.2(mat, Rowv=as.dendrogram(hc), symm=T, trace='none',
          col=rev(hmcol), margin=c(13,13))
dev.off()

