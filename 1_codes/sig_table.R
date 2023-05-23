rm(list = ls())
require(data.table)
library(tidyverse)

load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_fantom_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_ensembl_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_gencode_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_lncrnakb_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_noncode_gene.Rdata")

# ensembl <- rtracklayer::import("../annot/Homo_sapiens.GRCh38.104.gtf.gz")
# ensembl <- as.data.frame(ensembl)
# lncrna_ensembl_gene <- ensembl %>%
#   filter(type == "gene",gene_biotype == "lncRNA") %>%
#   select(start,end,gene_id,gene_name)%>%
#   distinct()
# save(lncrna_ensembl_gene,file = "lncrna_ensembl_gene.Rdata")
# 
# gencode <- rtracklayer::import("../annot/gencode.v38.annotation.gtf")
# gencode <- as.data.frame(gencode)
# lncrna_gencode_gene <- gencode %>% 
#   filter(type == c("gene"),gene_type == "lncRNA") %>%
#   select(start,end,gene_id,gene_name)
# save(lncrna_gencode_gene,file = "lncrna_gencode_gene.Rdata")
# 
# noncode <- read.table("../annot/NONCODEv6_hg38.lncAndGene.bed.gz")
# lncrna_noncode_gene <- noncode %>% select(start=V2,end=V3,gene_id=V4) %>% 
#   filter(str_detect(gene_id, 'NONHSAG')) %>%
#   distinct()
# save(lncrna_noncode_gene,file = "lncrna_noncode_gene.Rdata")
# 
# lncrnakb <- rtracklayer::import("../annot/lncRNAKB_hg38_lnconly_v7.gtf")
# lncrnakb <- as.data.frame(lncrnakb)
# lncrna_lncrnakb_gene <- lncrnakb %>% 
#   filter(type == "gene", GENE_TYPE == "lncRNA") %>%
#   select(start,end,gene_id,gene_name) %>% 
#   distinct()
# save(lncrna_lncrnakb_gene,file = "lncrna_lncrnakb_gene.Rdata")
# 
# fantom <- rtracklayer::import("../annot/FANTOM_CAT.lv3_robust.gtf.gz")
# fantom <- as.data.frame(fantom)
# lncrna_fantom_gene <- fantom %>% 
#   filter(type == "gene",geneSuperClass == "all_lncRNA") %>%
#   select(start,end,gene_id,gene_name) %>% 
#   distinct()
# save(lncrna_fantom_gene,file = "lncrna_fantom_gene.Rdata")

# ensembl <- rtracklayer::import("./annot/Homo_sapiens.GRCh38.104.gtf.gz")
# ensembl <- as.data.frame(ensembl)
# 
# mendelian_genes <- read.table("/rprojectnb2/pelosolab/yuxuan/glgc/mendelian_gene.txt")
# mendelian_genes <- mendelian_genes$V1 
# #mendelian_genes <- c("ABCA1", "ABCG5", "ANGPTL3", "APOA5", "APOB", 
# #                     "APOE", "CETP", "CYP27A1", "GPD1", "GPIHBP1", "LCAT", 
# #                     "LDLR", "LDLRAP1", "LIPA", "LIPC", "LMF1", "LPL", "MTTP", "PCSK9", "SAR1B", "SCARB1")
# mendelian_genes <- ensembl[ensembl$gene_name %in% mendelian_genes,]
# 
# mendelian_genes <- mendelian_genes %>%   
#   filter(type == "gene",gene_biotype == "protein_coding") %>% 
#   select(gene_name, seqnames, start,end)%>%
#   distinct() %>%
#   mutate(group = "mendelian gene")
# 
# 
# colnames(mendelian_genes)[c(1,2)] <- c("gene_id","chr")
# mendelian_genes$chr <- as.numeric(mendelian_genes$chr)

load("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/mendelian_genes.Rdata")

mendelian_genes$start_100k <- mendelian_genes$start - 100000
mendelian_genes$end_100k <- mendelian_genes$end + 100000

setDT(mendelian_genes)
setkey(mendelian_genes,chr,start_100k,end_100k)

path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results")
trait <- c("LDL_ADJ.norm","HDL.norm","TOTAL_ADJ.norm","TG_LOG.norm")
#trait <- "LDL_ADJ.norm"
annot <- c("gencode","fantom","lncrnakb","noncode")

for (i in trait){

  df <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/",i,"_",annot[1],".csv"))
  gwasResults <- inner_join(df, lncrna_gencode_gene)
  gwasResults_gencode <- gwasResults[,c(1,2,91,92,90)] 
  gwasResults_gencode$group <- "GENCODE"
  
  df <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/",i,"_",annot[2],".csv"))
  gwasResults <- inner_join(df, lncrna_fantom_gene)
  gwasResults_fantom <- gwasResults[,c(1,2,91,92,90)]
  gwasResults_fantom$group <- "FANTOM"
  
  df <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/",i,"_",annot[3],".csv"))
  gwasResults <- inner_join(df, lncrna_lncrnakb_gene)
  gwasResults_lncrnakb <- gwasResults[,c(1,2,91,92,90)]
  gwasResults_lncrnakb$group <- "LncrnaKB"
  
  df <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/",i,"_",annot[4],".csv"))
  gwasResults <- inner_join(df, lncrna_noncode_gene)
  gwasResults_noncode <- gwasResults[,c(1,2,91,92,90)]
  gwasResults_noncode$group <- "NONCODE"
  
  gwasResults <- rbind(gwasResults_gencode,gwasResults_fantom,gwasResults_lncrnakb,gwasResults_noncode)
  gwasResults$chr <- as.numeric(gsub("chr","",gwasResults$chr))
  gwasResults <- arrange(gwasResults, chr, start)
  sig_genes <- gwasResults[gwasResults$`STAAR-O` < 0.05/111550,]
  
  setDT(sig_genes)
  setkey(mendelian_genes,chr,start_100k,end_100k)
  setkey(sig_genes,chr,start,end)
  overlap <- foverlaps(sig_genes, mendelian_genes, type="any",nomatch = NA)
  write.csv(overlap,paste0("1_results/",i,"_sig_genes.csv"))
}


##############################



