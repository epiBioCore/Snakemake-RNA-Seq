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
KEGG_Module_species <- substr(species,1,3)

load_bioconductor_package(snakemake@input[["species_anno"]], species_pkg)


sig_genes <- read.csv(sig_genes)

enrich <- enrichMKEGG(gene = as.character(sig_genes$ENTREZID), 
          organism = KEGG_Module_species,
pAdjustMethod = "BH",
qvalueCutoff = 0.05)

if (!is.null(enrich)) {
KEGG_Module <- setReadable(enrich,OrgDb = get(species_pkg),keyType = "ENTREZID")

name <- paste(comp,"KEGG_Module",sep= "_")
assign(name,KEGG_Module)  
save(list=name,file = out)



if( nrow(KEGG_Module)>1) {
write.table(KEGG_Module,file=file.path(dir,paste0(comp,"_KEGG_Module.txt")), sep="\t", quote=F,row.names=F)

b <- barplot(KEGG_Module)
ggsave(b,file=file.path(dir,paste0(comp,"_KEGG_Module_barplot.png")),width = width,height = height,units = units)

b2 <- barplot(KEGG_Module,showCategory=12,drop=T)
ggsave(b2,file=file.path(dir,paste0(comp,"_KEGG_Module_barplot2.png")),width = width,height = height,units = units)

e <- emapplot(pairwise_termsim(KEGG_Module))
ggsave(e,file=file.path(dir,paste0(comp,"_KEGG_Module_emma.png")),width = width,height = height,units = units)

c <- cnetplot(KEGG_Module)
ggsave(c, file=file.path(dir,paste0(comp,"_KEGG_Module_cnet.png")),width = width,height = height,units = units)

d <- dotplot(KEGG_Module)
ggsave(d, file=file.path(dir,paste0(comp,"_KEGG_Module_dotplot.png")),width = width,height = height,units = units)


} else {
    text = paste("\n   No plot generated because no KEGG module was found significant\n",
         "       for comparison ",comp,"\n")
 p <- ggplot() + 
  annotate("text", x = 4, y = 25, size=4, label = text) + 
  theme_void()

 plots <- c("barplot","barplot2","emma","cnet","dotplot")  
 map(plots,~ggsave(p,file=file.path(dir,paste0(comp,"_KEGG_Module_",.x,".png"))))


}
} else {
    cat("No gene mapped to KEGG Module genes",file = file.path(dir,paste0(comp,"_KEGG_Module.txt")))
}