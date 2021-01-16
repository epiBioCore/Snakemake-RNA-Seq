log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(tidyverse)
library(ggrepel)
theme_set(theme_bw())
load(snakemake@input[["dds"]])

## function
plotPCA2 <- function (object, intgroup = "condition", ntop = 500, dim1= 1,dim2 = 2,returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }

  d <- data.frame( pca$x[, dim1],pca$x[, dim2], group = group, 
                  intgroup.df, name = colnames(object),descriptive_name = colData(object)$descriptive_name)
                
  dim_names <- paste0("PC",c(dim1,dim2))
  colnames(d)[1:2] <- dim_names
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(dim1,dim2)]
    return(d)
  }
  
  ggplot(data = d, aes_string(x = dim_names[1], y = dim_names[2], color = "group")) + 
    geom_point(size = 3) + xlab(paste0(dim_names[1],": ", round(percentVar[dim1] * 
                                                        100), "% variance")) + ylab(paste0(dim_names[2],": ", round(percentVar[dim2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}





rld <- rlog(dds)
save(rld,file=snakemake@output[["rld"]])

dims <- list(c(1,2),c(2,3),c(1,3))

colData(rld)
size = as.integer(snakemake@params[["lab_size"]])
size
lapply(dims,function(d) {
    dim1<-paste0("PC",d[1])
    dim2<-paste0("PC",d[2])

    lapply(snakemake@params[["col"]],function(v) {
        p <- plotPCA2(rld,intgroup = v,dim1 = d[1],dim2=d[2])

       if (snakemake@params[["labs"]]) {
          p <-  p +
           geom_text_repel(aes(label = descriptive_name),size = as.integer(snakemake@params[["lab_size"]]))
       }
       ggsave(p,file=file.path(dirname(snakemake@input[["dds"]]),paste0("PCA_",dim1,"vs",dim2,"_colored_by_",v,".pdf")))         
    })
})   