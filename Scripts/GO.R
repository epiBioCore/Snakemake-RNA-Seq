log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggnewscale)
source(file.path(snakemake@scriptdir, 'Common.R') )
species_pkg <- snakemake@params[["species"]]


load_bioconductor_package(snakemake@input[["species_anno"]], species_pkg)

ont= snakemake@wildcards$GO

sig_genes <- snakemake@input[["sig"]]
genelist <- sub("SIG","ALL",sig_genes)
out <- snakemake@output[["object"]]
dir <- snakemake@params[["out"]]
comp= snakemake@wildcards$contrast
height <- snakemake@params[["h"]]
width <- snakemake@params[["w"]]
units <- snakemake@params[["u"]]

sig_genes <- read.csv(sig_genes)
genelist <- read.csv(genelist)

enrich <- enrichGO(gene =as.character(sig_genes$ENTREZID),
         universe =as.character(genelist$ENTREZID),
         OrgDb=get(species_pkg), ont=gsub("GO_","",ont),
         pAdjustMethod = "BH",
         pvalueCutoff  = 0.01,
         qvalueCutoff  = 0.05, readable =TRUE)


enrich
name <- paste(comp,ont,sep= "_")
assign(name,enrich)  
save(list=name,file = out)

save(enrich,file=out)
if(nrow(enrich)>1) {
write.table(enrich,file=file.path(dir,paste0(comp,"_",ont,".txt")), sep="\t", quote=F,row.names=F)


b <- barplot(enrich)
ggsave(b,file=file.path(dir,paste0(comp,"_",ont,"_barplot.png")),width = width,height = height,units = units)

b2 <- barplot(enrich,showCategory=12,drop=T)
ggsave(b2,file=file.path(dir,paste0(comp,"_",ont,"_barplot2.png")),width = width,height = height,units = units)


e <- emapplot(pairwise_termsim(enrich))
ggsave(e,file=file.path(dir,paste0(comp,"_",ont,"_emma.png")),width = width,height = height,units = units)

c <- cnetplot(enrich)
ggsave(c,file=file.path(dir,paste0(comp,"_",ont,"_cnet.png")),width = width,height = height,units = units)

d <- dotplot(enrich)
ggsave(d,file=file.path(dir,paste0(comp,"_",ont,"_dotplot.png")),width = width,height = height,units = units)


} else {
    text = paste("\n   No plot generated because no", ont, " terms found significant \n",
         "for comparison ",comp,"\n")
         
 p <- ggplot() + 
  annotate("text", x = 4, y = 25, size=4, label = text) + 
  theme_void()

 plots <- c("barplot","barplot2","emma","cnet","dotplot")  
 map(plots,~ggsave(p,file=file.path(dir,paste0(comp,"_",ont,"_",.x,".png"))))


}