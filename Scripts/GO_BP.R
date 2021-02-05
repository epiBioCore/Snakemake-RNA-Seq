log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

source(file.path(snakemake@scriptdir, 'Common.R') )
species_pkg <- snakemake@params[["species"]]


load_bioconductor_package(snakemake@input[["species_anno"]], species_pkg)

ont= snakemake@wildcards$GO
sig_genes <- snakemake@input[["sig"]]
genelist <- sub("SIG","ALL",sig_genes)
out <- snakemake@output[["object"]]
dir <- snakemake@params[["out"]]
comp= snakemake@wildcards$contrast
height <- snakemake@params[["height"]]
width <- snakemake@params[["width"]]
units <- snakemake@params[["units"]]


enrich <- enrichGO(gene =sig_genes$ENTREZID,
         universe =genelists$ENTREZID,
         OrgDb=get(species_pkg), ont=ont,
         pAdjustMethod = "BH",
         pvalueCutoff  = 0.01,
         qvalueCutoff  = 0.05, readable =TRUE)


save(enrich,file=object_out)
if(nrow(enrich)>0) {
write.table(enrich,file=file.path(out), sep="\t", quote=F,row.names=F)

barplot(enrich)
ggsave(file=file.path(out,paste0(comp,"_",ont,"_barplot.png"),width = width,height = height,units = units)

barplot(enrich,showCategory=12,drop=T)
ggsave(file=file.path(out,paste0(comp,"_",ont,"_barplot2.png"),width = width,height = height,units = units)

emapplot(enrich)
ggsave(file=file.path(out,paste0(comp,"_",ont,"_emma.png"),width = width,height = height,units = units)

cnetplot(enrich)
ggsave(file=file.path(out,paste0(comp,"_",ont,"_cnet.png"),width = width,height = height,units = units)

dotplot(enrich)
ggsave(file=file.path(out,paste0(comp,"_",ont,"_dotplot.png"),width = width,height = height,units = units)


} else {
    text = paste("\n   No plot generated, because no", GO, " terms found significant.\n",
         "       for comparison ",comp,"\n")
 p <- ggplot() + 
  annotate("text", x = 4, y = 25, size=8, label = text) + 
  theme_void()

 plots <- c("barplot","barplot2","emma","cnet","dotplot")  
 map(plots,~ggsave(p,file=file.path(out,paste0(comp,"_",ont,"_",.x,".png"),width = width,height = height,units = units))


}