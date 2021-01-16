out1 <- 
out2 <- 
out3
out4
out5
 
in <- 
load(in)
barplot(enrich,showCategory=12)
ggsave(out1)

barplot(enrich,showCategory=12,drop=T))
ggsave(out2)

emapplot(enrich)
ggsave(out3)

cnetplot(enrich)
ggsave(out4)

dotplot(enrich)
ggsave(out5)