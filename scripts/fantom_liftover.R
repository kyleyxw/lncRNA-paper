library(dplyr)
library(stringr)
library(GenomicRanges)
library(rtracklayer)
# #####################################################################


## Fantom gene
fantom <- rtracklayer::import("FANTOM_CAT.lv3_robust.gtf.gz")
fantom <- as.data.frame(fantom)
lncrna_fantom_gene <- fantom %>% 
  filter(type == "gene",geneSuperClass == "all_lncRNA") %>%
  select(seqnames,start,end,gene_id,gene_name) %>% 
  distinct()

write.table(lncrna_fantom_gene , "lncrna_fantom_gene.bed", row.names = F,col.names = F,quote = F,sep = '\t')
## lift over from ucsc web

lncrna_fantom_gene_hg38 <- read.table(file = "lncrna_fantom_gene_hg38.bed",stringsAsFactors = F)
colnames(lncrna_fantom_gene_hg38) <- c(colnames(lncrna_fantom_gene))

#lncrna_fantom_gene_hg38 <- lncrna_fantom_gene_hg38[!lncrna_fantom_gene_hg38$gene_name %in% lncrna_gencode_gene$gene_name,]
lncrna_fantom_gene <- lncrna_fantom_gene_hg38
lncrna_fantom_gene$gene_type <- "lncRNA"
save(lncrna_fantom_gene,file="lncrna_fantom_gene.Rdata")

chr <- c(paste0("chr",1:22))
for (i in chr){
  lncrna_fantom_gene <- lncrna_fantom_gene_hg38 %>% 
    filter(seqnames == i)
  save(lncrna_fantom_gene,file=paste0("lncrna_fantom_gene_",i,".Rdata"))
}


