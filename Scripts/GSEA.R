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

genelist <- snakemake@input[["de"]]
out <- snakemake@output[["object"]]
dir <- snakemake@params[["out"]]
comp= snakemake@wildcards$contrast
height <- snakemake@params[["h"]]
width <- snakemake@params[["w"]]
units <- snakemake@params[["u"]]


genelist <- read.csv(genelist)
########################################
  
top<-genelist[!is.na(genelist$ENTREZID),]
  
top<-top[!duplicated(top$ENTREZID),]
  
GSEA_list=top$log2FoldChange
  
names(GSEA_list)<-as.vector(as.character(top$ENTREZID))
  
GSEA_list=sort(GSEA_list, decreasing=T)
  
gse <- gseGO(geneList     = GSEA_list,
        
        OrgDb        = get(species_pkg),
        
        ont          = "ALL",
        
        nPerm        = 1000,
        
        minGSSize    = 100,
        
        maxGSSize    = 500,
        
        pvalueCutoff = 0.05,

        
        verbose      = FALSE)

name <- paste(comp,"GSEA",sep= "_")
assign(name,gse)  
save(list=name,file = out)
gse <- setReadable(gse,get(species_pkg),keyType = "ENTREZID")


if(nrow(gse)>1) {
write.table(gse,file=file.path(dir,paste0(comp,"GSEA.txt")), sep="\t", quote=F,row.names=F)


e <- emapplot(pairwise_termsim(gse))
ggsave(e,file=file.path(dir,paste0(comp,"_GSEA_emma.png")),width = width,height = height,units = units)

c <- cnetplot(gse)
ggsave(c, file=file.path(dir,paste0(comp,"_GSEA_cnet.png")),width = width,height = height,units = units)

d <- dotplot(gse)
ggsave(d,file=file.path(dir,paste0(comp,"_GSEA_dotplot.png")),width = width,height = height,units = units)


} else {
    text = paste("\n   No plot generated because no GO terms found significant\n",
         "       for comparison ",comp,"\n")
 p <- ggplot() + 
  annotate("text", x = 4, y = 25, size=4, label = text) + 
  theme_void()

 plots <- c("emma","cnet","dotplot")  
 map(plots,~ggsave(p,file=file.path(dir,paste0(comp,"_gsea_",.x,"_plot.png"))))


}

