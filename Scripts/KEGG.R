log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggnewscale)
source(file.path(snakemake@scriptdir, 'Common.R') )





sig_genes <- snakemake@input[["sig"]]
out <- snakemake@output[["object"]]
dir <- snakemake@params[["out"]]
comp= snakemake@wildcards$contrast
height <- snakemake@params[["h"]]
width <- snakemake@params[["w"]]
units <- snakemake@params[["u"]]
species_pkg <- snakemake@params[["species_anno"]]
species <- snakemake@params[["species"]]
kegg_species <- substr(species,1,3)


load_bioconductor_package(snakemake@input[["species_anno"]], species_pkg)


sig_genes <- read.csv(sig_genes)

enrich <- enrichKEGG(gene = as.character(sig_genes$ENTREZID), 
          organism = kegg_species,
pAdjustMethod = "BH",
qvalueCutoff = 0.05)

KEGG <- setReadable(enrich,OrgDb = get(species_pkg),keyType = "ENTREZID")

name <- paste(comp,"KEGG",sep= "_")
assign(name,KEGG)  
save(list=name,file = out)



if(nrow(KEGG)>1) {
write.table(KEGG,file=file.path(dir,paste0(comp,"_KEGG.txt")), sep="\t", quote=F,row.names=F)

b <- barplot(KEGG)
ggsave(b,file=file.path(dir,paste0(comp,"_KEGG_barplot.png")),width = width,height = height,units = units)

b2 <- barplot(KEGG,showCategory=12,drop=T)
ggsave(b2,file=file.path(dir,paste0(comp,"_KEGG_barplot2.png")),width = width,height = height,units = units)

e <- emapplot(pairwise_termsim(KEGG))
ggsave(e,file=file.path(dir,paste0(comp,"_KEGG_emma.png")),width = width,height = height,units = units)

c <- cnetplot(KEGG)
ggsave(c, file=file.path(dir,paste0(comp,"_KEGG_cnet.png")),width = width,height = height,units = units)

d <- dotplot(KEGG)
ggsave(d, file=file.path(dir,paste0(comp,"_KEGG_dotplot.png")),width = width,height = height,units = units)


} else {
    text = paste("\n   No plot generated, because no terms found significant.\n",
         "       for comparison ",comp,"\n")
 p <- ggplot() + 
  annotate("text", x = 4, y = 25, size=4, label = text) + 
  theme_void()

 plots <- c("barplot","barplot2","emma","cnet","dotplot")  
 map(plots,~ggsave(p,file=file.path(dir,paste0(comp,"_KEGG_",.x,".png"))))


}