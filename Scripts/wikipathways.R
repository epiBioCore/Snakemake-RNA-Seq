log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggnewscale)
library(rWikiPathways)
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

l <- listOrganisms()

### grepl(substr(species,1,1,),species) & grepl(sub(".","",species))
### input species is in format hsapiens
### Species in wikipathways are in format: Homo sapiens
## so extracting and capitalizing first letter of input species: h -> H
## and also the species name (sapiens)
## grepl is finding the wiki organism that starts with the captilized first letter and species name
first <- paste0("^",toupper(substr(species,1,1)))

organism <- l[grepl(first,l) & grepl(sub(".","",species),l)]


load_bioconductor_package(snakemake@input[["species_anno"]], species_pkg)


sig_genes <- read.csv(sig_genes)
##

wp.gmt <- rWikiPathways::downloadPathwayArchive(organism=organism, format = "gmt")
wp.gmt
wp2gene <- clusterProfiler::read.gmt(wp.gmt)
wp2gene
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME

ewp <- enricher(sig_genes$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name,pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)

ewp <- setReadable(ewp,get(species_pkg), keyType = "ENTREZID")

name <- paste(comp,"wikipathways",sep= "_")
assign(name,ewp)  
save(list=name,file = out)



if(nrow(ewp)>1) {
write.table(ewp,file=file.path(dir,paste0(comp,"_wikipathways.txt")), sep="\t", quote=F,row.names=F)

b <- barplot(ewp)
ggsave(b,file=file.path(dir,paste0(comp,"_wikipathways_barplot.png")),width = width,height = height,units = units)

b2 <- barplot(ewp,showCategory=12,drop=T)
ggsave(b2,file=file.path(dir,paste0(comp,"_wikipathways_barplot2.png")),width = width,height = height,units = units)

e <- emapplot(pairwise_termsim(ewp))
ggsave(e,file=file.path(dir,paste0(comp,"_wikipathways_emma.png")),width = width,height = height,units = units)

c <- cnetplot(ewp)
ggsave(c, file=file.path(dir,paste0(comp,"_ewp_cnet.png")),width = width,height = height,units = units)

d <- dotplot(ewp)
ggsave(d, file=file.path(dir,paste0(comp,"_ewp_dotplot.png")),width = width,height = height,units = units)


} else {
    text = paste("\n   No plot generated because no wikipathways were significant\n",
         "       for comparison ",comp,"\n")
 p <- ggplot() + 
  annotate("text", x = 4, y = 25, size=4, label = text) + 
  theme_void()

 plots <- c("barplot","barplot2","emma","cnet","dotplot")  
 map(plots,~ggsave(p,file=file.path(dir,paste0(comp,"_wikipathways_",.x,".png")))) 


}
