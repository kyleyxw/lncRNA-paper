rm(list = ls())
library(meta)
library(tidyverse)
library(qqman)
library(CMplot)

#####################
fantom <- rtracklayer::import("../annot/FANTOM_CAT.lv3_robust.gtf.gz")
fantom <- as.data.frame(fantom)
lncrna_fantom_gene <- fantom %>% 
  filter(type == "gene",geneSuperClass == "all_lncRNA") %>%
  select(start,end,gene_id,gene_name) %>% 
  distinct()

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/LDL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_fantom_gene_",full.names = TRUE)
ldl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
ldl <- do.call(rbind, ldl)
ldl <- inner_join(ldl, lncrna_fantom_gene) 
ldl$chr <- as.numeric(str_remove(ldl$chr, "chr"))
ldl$`STAAR-O` <- as.numeric(as.character(ldl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/HDL.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_fantom_gene_",full.names = TRUE)
hdl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
hdl <- do.call(rbind, hdl)
hdl <- inner_join(hdl, lncrna_fantom_gene)
hdl$chr <- as.numeric(str_remove(hdl$chr, "chr"))
hdl$`STAAR-O` <- as.numeric(as.character(hdl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TOTAL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_fantom_gene_",full.names = TRUE)
total <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
total <- do.call(rbind, total)
total <- inner_join(total, lncrna_fantom_gene) 
total$chr <- as.numeric(str_remove(total$chr, "chr"))
total$`STAAR-O` <- as.numeric(as.character(total$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TG_LOG.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_fantom_gene_",full.names = TRUE)
tg <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
tg <- do.call(rbind, tg)
tg <- inner_join(tg, lncrna_fantom_gene) 
tg$chr <- as.numeric(str_remove(tg$chr, "chr"))
tg$`STAAR-O` <- as.numeric(as.character(tg$`STAAR-O`))

df <- ldl[,c(1,2,91,90)] %>%
  inner_join(hdl[,c(1,2,91,90)],by = c("chr", "start", "gene_id"), suffix = c("_LDL", "_HDL")) %>%
  inner_join(total[,c(1,2,91,90)],by = c("chr", "start", "gene_id")) %>%
  inner_join(tg[,c(1,2,91,90)],by = c("chr", "start","gene_id"), suffix = c("_Total", "_TG"))

SNPs <- list(
  df$gene_id[df$`STAAR-O_LDL` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_HDL` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_Total` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_TG` < 0.05 / nrow(df)]
)


CMplot(df,type="p", plot.type="b",multracks=TRUE,threshold=0.05/nrow(df),threshold.lty=2,threshold.col=c("grey"), 
       highlight.cex = 1, highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4, highlight.pch=16,
       amplify=FALSE, bin.size=1e6,signal.col=NULL, file="jpg",memo="fantom",
       dpi=300,file.output=TRUE,verbose=TRUE)


#####################
gencode <- rtracklayer::import("../annot/gencode.v38.annotation.gtf")
gencode <- as.data.frame(gencode)

lncrna_gencode_gene <- gencode %>% 
  filter(type == c("gene"),gene_type == "lncRNA") %>%
  select(start,end,gene_id,gene_name)

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/LDL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_gencode_gene_",full.names = TRUE)
ldl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
ldl <- do.call(rbind, ldl)
ldl <- inner_join(ldl, lncrna_gencode_gene) 
ldl$chr <- as.numeric(str_remove(ldl$chr, "chr"))
ldl$`STAAR-O` <- as.numeric(as.character(ldl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/HDL.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_gencode_gene_",full.names = TRUE)
hdl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
hdl <- do.call(rbind, hdl)
hdl <- inner_join(hdl, lncrna_gencode_gene)
hdl$chr <- as.numeric(str_remove(hdl$chr, "chr"))
hdl$`STAAR-O` <- as.numeric(as.character(hdl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TOTAL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_gencode_gene_",full.names = TRUE)
total <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
total <- do.call(rbind, total)
total <- inner_join(total, lncrna_gencode_gene) 
total$chr <- as.numeric(str_remove(total$chr, "chr"))
total$`STAAR-O` <- as.numeric(as.character(total$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TG_LOG.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_gencode_gene_",full.names = TRUE)
tg <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
tg <- do.call(rbind, tg)
tg <- inner_join(tg, lncrna_gencode_gene) 
tg$chr <- as.numeric(str_remove(tg$chr, "chr"))
tg$`STAAR-O` <- as.numeric(as.character(tg$`STAAR-O`))

df <- ldl[,c(93,2,91,90)] %>%
  inner_join(hdl[,c(93,2,91,90)],by = c("chr", "start", "gene_name"), suffix = c("_LDL", "_HDL")) %>%
  inner_join(total[,c(93,2,91,90)],by = c("chr", "start", "gene_name")) %>%
  inner_join(tg[,c(93,2,91,90)],by = c("chr", "start","gene_name"), suffix = c("_Total", "_TG"))

SNPs <- list(
  df$gene_name[df$`STAAR-O_LDL` < 0.05 / nrow(df)],
  df$gene_name[df$`STAAR-O_HDL` < 0.05 / nrow(df)],
  df$gene_name[df$`STAAR-O_Total` < 0.05 / nrow(df)],
  df$gene_name[df$`STAAR-O_TG` < 0.05 / nrow(df)]
)


CMplot(df,type="p", plot.type="b",multracks=TRUE,threshold=0.05/nrow(df),threshold.lty=2,threshold.col=c("grey"), 
       highlight.cex = 1, highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4, highlight.pch=16,
       amplify=FALSE, bin.size=1e6,signal.col=NULL, file="jpg",memo="gencode",
       dpi=300,file.output=TRUE,verbose=TRUE)

#####################
ensembl <- rtracklayer::import("../annot/Homo_sapiens.GRCh38.104.gtf.gz")
ensembl <- as.data.frame(ensembl)

lncrna_ensembl_gene <- ensembl %>%
  filter(type == "gene",gene_biotype == "lncRNA") %>%
  select(start,end,gene_id,gene_name)%>%
  distinct()

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/LDL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_ensembl_gene_",full.names = TRUE)
ldl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
ldl <- do.call(rbind, ldl)
ldl <- inner_join(ldl, lncrna_ensembl_gene) 
ldl$chr <- as.numeric(str_remove(ldl$chr, "chr"))
ldl$`STAAR-O` <- as.numeric(as.character(ldl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/HDL.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_ensembl_gene_",full.names = TRUE)
hdl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
hdl <- do.call(rbind, hdl)
hdl <- inner_join(hdl, lncrna_ensembl_gene)
hdl$chr <- as.numeric(str_remove(hdl$chr, "chr"))
hdl$`STAAR-O` <- as.numeric(as.character(hdl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TOTAL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_ensembl_gene_",full.names = TRUE)
total <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
total <- do.call(rbind, total)
total <- inner_join(total, lncrna_ensembl_gene) 
total$chr <- as.numeric(str_remove(total$chr, "chr"))
total$`STAAR-O` <- as.numeric(as.character(total$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TG_LOG.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_ensembl_gene_",full.names = TRUE)
tg <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
tg <- do.call(rbind, tg)
tg <- inner_join(tg, lncrna_ensembl_gene) 
tg$chr <- as.numeric(str_remove(tg$chr, "chr"))
tg$`STAAR-O` <- as.numeric(as.character(tg$`STAAR-O`))

df <- ldl[,c(1,2,91,90)] %>%
  inner_join(hdl[,c(1,2,91,90)],by = c("chr", "start", "gene_id"), suffix = c("_LDL", "_HDL")) %>%
  inner_join(total[,c(1,2,91,90)],by = c("chr", "start", "gene_id")) %>%
  inner_join(tg[,c(1,2,91,90)],by = c("chr", "start","gene_id"), suffix = c("_Total", "_TG"))

SNPs <- list(
  df$gene_id[df$`STAAR-O_LDL` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_HDL` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_Total` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_TG` < 0.05 / nrow(df)]
)


CMplot(df,type="p", plot.type="b",multracks=TRUE,threshold=0.05/nrow(df),threshold.lty=2,threshold.col=c("grey"), 
       highlight.cex = 1, highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4, highlight.pch=16,
       amplify=FALSE, bin.size=1e6,signal.col=NULL, file="jpg",memo="ensembl",
       dpi=300,file.output=TRUE,verbose=TRUE)

#####################
noncode <- read.table("../annot/NONCODEv6_hg38.lncAndGene.bed.gz")

lncrna_noncode_gene <- noncode %>% select(start=V2,end=V3,gene_id=V4) %>% 
  filter(str_detect(gene_id, 'NONHSAG')) %>%
  distinct()

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/LDL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_noncode_gene_",full.names = TRUE)
ldl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
ldl <- do.call(rbind, ldl)
ldl <- inner_join(ldl, lncrna_noncode_gene) 
ldl$chr <- as.numeric(str_remove(ldl$chr, "chr"))
ldl$`STAAR-O` <- as.numeric(as.character(ldl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/HDL.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_noncode_gene_",full.names = TRUE)
hdl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
hdl <- do.call(rbind, hdl)
hdl <- inner_join(hdl, lncrna_noncode_gene)
hdl$chr <- as.numeric(str_remove(hdl$chr, "chr"))
hdl$`STAAR-O` <- as.numeric(as.character(hdl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TOTAL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_noncode_gene_",full.names = TRUE)
total <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
total <- do.call(rbind, total)
total <- inner_join(total, lncrna_noncode_gene) 
total$chr <- as.numeric(str_remove(total$chr, "chr"))
total$`STAAR-O` <- as.numeric(as.character(total$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TG_LOG.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_noncode_gene_",full.names = TRUE)
tg <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
tg <- do.call(rbind, tg)
tg <- inner_join(tg, lncrna_noncode_gene) 
tg$chr <- as.numeric(str_remove(tg$chr, "chr"))
tg$`STAAR-O` <- as.numeric(as.character(tg$`STAAR-O`))

df <- ldl[,c(1,2,91,90)] %>%
  full_join(hdl[,c(1,2,91,90)],by = c("chr", "start", "gene_id"), suffix = c("_LDL", "_HDL")) %>%
  full_join(total[,c(1,2,91,90)],by = c("chr", "start", "gene_id")) %>%
  full_join(tg[,c(1,2,91,90)],by = c("chr", "start","gene_id"), suffix = c("_Total", "_TG"))

df[is.na(df)] <- 1

SNPs <- list(
  df$gene_id[df$`STAAR-O_LDL` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_HDL` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_Total` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_TG` < 0.05 / nrow(df)]
)


CMplot(df,type="p", plot.type="b",multracks=TRUE,threshold=0.05/nrow(df),threshold.lty=2,threshold.col=c("grey"), 
       highlight.cex = 1, highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4, highlight.pch=16,
       amplify=FALSE, bin.size=1e6,signal.col=NULL, file="jpg",memo="noncode",
       dpi=300,file.output=TRUE,verbose=TRUE)

#####################
lncrnakb <- rtracklayer::import("../annot/lncRNAKB_hg38_lnconly_v7.gtf")
lncrnakb <- as.data.frame(lncrnakb)

lncrna_lncrnakb_gene <- lncrnakb %>% 
  filter(type == "gene", GENE_TYPE == "lncRNA") %>%
  select(start,end,gene_id,gene_name) %>% 
  distinct()

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/LDL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_lncrnakb_gene_",full.names = TRUE)
ldl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
ldl <- do.call(rbind, ldl)
ldl <- inner_join(ldl, lncrna_lncrnakb_gene) 
ldl$chr <- as.numeric(str_remove(ldl$chr, "chr"))
ldl$`STAAR-O` <- as.numeric(as.character(ldl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/HDL.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_lncrnakb_gene_",full.names = TRUE)
hdl <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
hdl <- do.call(rbind, hdl)
hdl <- inner_join(hdl, lncrna_lncrnakb_gene)
hdl$chr <- as.numeric(str_remove(hdl$chr, "chr"))
hdl$`STAAR-O` <- as.numeric(as.character(hdl$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TOTAL_ADJ.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_lncrnakb_gene_",full.names = TRUE)
total <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
total <- do.call(rbind, total)
total <- inner_join(total, lncrna_lncrnakb_gene) 
total$chr <- as.numeric(str_remove(total$chr, "chr"))
total$`STAAR-O` <- as.numeric(as.character(total$`STAAR-O`))

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/TG_LOG.norm/"
filelist <- list.files(path,pattern="topmed_lncRNA_lncrnakb_gene_",full.names = TRUE)
tg <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
tg <- do.call(rbind, tg)
tg <- inner_join(tg, lncrna_lncrnakb_gene) 
tg$chr <- as.numeric(str_remove(tg$chr, "chr"))
tg$`STAAR-O` <- as.numeric(as.character(tg$`STAAR-O`))

df <- ldl[,c(1,2,91,90)] %>%
  full_join(hdl[,c(1,2,91,90)],by = c("chr", "start", "gene_id"), suffix = c("_LDL", "_HDL")) %>%
  full_join(total[,c(1,2,91,90)],by = c("chr", "start", "gene_id")) %>%
  full_join(tg[,c(1,2,91,90)],by = c("chr", "start","gene_id"), suffix = c("_Total", "_TG"))

df[is.na(df)] <- 1

SNPs <- list(
  df$gene_id[df$`STAAR-O_LDL` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_HDL` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_Total` < 0.05 / nrow(df)],
  df$gene_id[df$`STAAR-O_TG` < 0.05 / nrow(df)]
)


CMplot(df,type="p", plot.type="b",multracks=TRUE,threshold=0.05/nrow(df),threshold.lty=2,threshold.col=c("grey"), 
       highlight.cex = 1, highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4, highlight.pch=16,
       amplify=FALSE, bin.size=1e6,signal.col=NULL, file="jpg",memo="lncrnakb",
       dpi=300,file.output=TRUE,verbose=TRUE)