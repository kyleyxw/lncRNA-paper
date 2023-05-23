rm(list = ls())
require(data.table)
library(tidyverse)

paste2 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

# load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_fantom_gene.Rdata")
# load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_ensembl_gene.Rdata")
# load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_gencode_gene.Rdata")
# load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_lncrnakb_gene.Rdata")
# load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_noncode_gene.Rdata")

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
# 
# save(mendelian_genes, file = "/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/mendelian_genes.Rdata")
setwd("/rprojectnb2/pelosolab/yuxuan/lncRNA/")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/mendelian_genes.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/glgc/trans_ancestry_assign_ind.Rdata")

mendelian_genes$start_100k <- mendelian_genes$start - 100000
mendelian_genes$end_100k <- mendelian_genes$end + 100000

setDT(mendelian_genes)
setkey(mendelian_genes,chr,start_100k,end_100k)

path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results")
filelist <- list.files(path,full.names = TRUE)


#trait <- c("LDL_ADJ.norm","HDL.norm","TOTAL_ADJ.norm","TG_LOG.norm")
trait <- "LDL_ADJ.norm_cond"
annot <- c("gencode","fantom","lncrnakb","noncode")


gencode_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[1],".csv"))
fantom_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[2],".csv"))
#fantom_lncrna <- fantom_lncrna[-grep("ENSG", fantom_lncrna$gene_id),]
lncrnakb_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[3],".csv"))
noncode_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[4],".csv"))

lncrna_genes <- rbind(gencode_lncrna,fantom_lncrna,lncrnakb_lncrna,noncode_lncrna)
lncrna_genes <- lncrna_genes[,c("gene_id","STAAR-O")]

ldl_sig_genes <- read_csv("1_results/LDL_ADJ.norm_sig_genes.csv")
ldl_sig_genes <- ldl_sig_genes[,c(9,13,2,10:11,3,12)]
colnames(ldl_sig_genes) <- c("gene_id","Annot","chr","start","end","Mendelian","STAAR-O")

ldl_sig_genes <- left_join(ldl_sig_genes,lncrna_genes,by = c("gene_id"))

ldl_sig_genes_cond <- rbind(ldl_sig_genes[is.na(ldl_sig_genes$Mendelian),],
                            aggregate(Mendelian ~ gene_id + Annot + chr + start + end + `STAAR-O.x` + `STAAR-O.y`,
                                      data = ldl_sig_genes, function(x) paste2(x, collapse=";")))
ldl_sig_genes_cond <- arrange(ldl_sig_genes_cond,chr,start)
colnames(ldl_sig_genes_cond)[c(7,8)] <- c("P","P_cond")

genes_bed <- ldl_sig_genes_cond[,c("chr","start", "end", "gene_id")]
genes_bed$chr <- paste0("chr",genes_bed$chr)
write.table(genes_bed, "2_results/genes_all.bed", row.names = F,col.names = F,quote = F,sep = '\t')

load("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/lipid_known_loci_info_rsID.Rdata")
snps_bed <- known_loci_info[,c("CHR","POS","POS","rsID")]
snps_bed <- snps_bed[!duplicated(snps_bed), ]
snps_bed$CHR <- paste0("chr",snps_bed$CHR)
write.table(snps_bed , "2_results/snps_all.bed", row.names = F,col.names = F,quote = F,sep = '\t')

setwd("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results")
Sys.setenv(PATH="/usr/bin/")
system("sort -k1,1 -k2,2n genes_all.bed > genes_all_sort.bed")
system("sort -k1,1 -k2,2n snps_all.bed > snps_all_sort.bed")

Sys.setenv(PATH="/share/pkg.7/bedtools/2.30.0/install/bin/")
system("bedtools closest -a genes_all_sort.bed -b snps_all_sort.bed  -t all -d -k 1 > snps_genes.bed")
Sys.setenv(PATH="/usr/bin/")

snps_genes <- read.table(file = "snps_genes.bed")
snps_genes <- snps_genes[,c("V4","V8","V7","V9")]
colnames(snps_genes) <- c("gene_id","rsID","pos","D")
ldl_sig_genes_cond <- left_join(ldl_sig_genes_cond,snps_genes) 

ldl_sig_genes_cond <- left_join(ldl_sig_genes_cond,
                                trans_ancestry_assign_ind[,c(1:2)],
                                by=c("rsID"="rsid_dbSNP150"))

ldl_sig_genes_cond <- ldl_sig_genes_cond %>%
  group_by(gene_id,Annot,chr,start,end,Mendelian,P,P_cond,rsID,pos,D) %>% summarise(trait = paste2(trait, collapse=","))
colnames(ldl_sig_genes_cond)[12] <- "glgc_trait"

ldl_sig_genes_cond <- ldl_sig_genes_cond[ldl_sig_genes_cond$P_cond < 0.05/length(unique(ldl_sig_genes_cond$gene_id)),]
ldl_sig_genes_cond <- arrange(ldl_sig_genes_cond,chr, start)

write_csv(ldl_sig_genes_cond,file = "/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/LDL_ADJ.norm_sig_genes_cond.csv")


##############################
setwd("/rprojectnb2/pelosolab/yuxuan/lncRNA/")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/mendelian_genes.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/glgc/trans_ancestry_assign_ind.Rdata")

mendelian_genes$start_100k <- mendelian_genes$start - 100000
mendelian_genes$end_100k <- mendelian_genes$end + 100000

setDT(mendelian_genes)
setkey(mendelian_genes,chr,start_100k,end_100k)

path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results")
trait <- "HDL.norm_cond"
annot <- c("gencode","fantom","lncrnakb","noncode")


gencode_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[1],".csv"))
fantom_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[2],".csv"))
#fantom_lncrna <- fantom_lncrna[-grep("ENSG", fantom_lncrna$gene_id),]
lncrnakb_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[3],".csv"))
noncode_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[4],".csv"))

lncrna_genes <- rbind(gencode_lncrna,fantom_lncrna,lncrnakb_lncrna,noncode_lncrna)
lncrna_genes <- lncrna_genes[,c("gene_id","STAAR-O")]

hdl_sig_genes <- read_csv("1_results/HDL.norm_sig_genes.csv")
hdl_sig_genes <- hdl_sig_genes[,c(9,13,2,10:11,3,12)]
colnames(hdl_sig_genes) <- c("gene_id","Annot","chr","start","end","Mendelian","STAAR-O")

hdl_sig_genes <- left_join(hdl_sig_genes,lncrna_genes,by = c("gene_id"))

hdl_sig_genes_cond <- rbind(hdl_sig_genes[is.na(hdl_sig_genes$Mendelian),],
                            aggregate(Mendelian ~ gene_id + Annot + chr + start + end + `STAAR-O.x` + `STAAR-O.y`,
                                      data = hdl_sig_genes, function(x) paste2(x, collapse=";")))
hdl_sig_genes_cond <- arrange(hdl_sig_genes_cond,chr,start)
colnames(hdl_sig_genes_cond)[c(7,8)] <- c("P","P_cond")

genes_bed <- hdl_sig_genes_cond[,c("chr","start", "end", "gene_id")]
genes_bed$chr <- paste0("chr",genes_bed$chr)
write.table(genes_bed, "2_results/genes_all.bed", row.names = F,col.names = F,quote = F,sep = '\t')

load("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/lipid_known_loci_info_rsID.Rdata")
snps_bed <- known_loci_info[,c("CHR","POS","POS","rsID")]
snps_bed <- snps_bed[!duplicated(snps_bed), ]
snps_bed$CHR <- paste0("chr",snps_bed$CHR)
write.table(snps_bed , "2_results/snps_all.bed", row.names = F,col.names = F,quote = F,sep = '\t')

setwd("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results")
Sys.setenv(PATH="/usr/bin/")
system("sort -k1,1 -k2,2n genes_all.bed > genes_all_sort.bed")
system("sort -k1,1 -k2,2n snps_all.bed > snps_all_sort.bed")

Sys.setenv(PATH="/share/pkg.7/bedtools/2.30.0/install/bin/")
system("bedtools closest -a genes_all_sort.bed -b snps_all_sort.bed  -t all -d -k 1 > snps_genes.bed")
Sys.setenv(PATH="/usr/bin/")

snps_genes <- read.table(file = "snps_genes.bed")
snps_genes <- snps_genes[,c("V4","V8","V7","V9")]
colnames(snps_genes) <- c("gene_id","rsID","pos","D")
hdl_sig_genes_cond <- left_join(hdl_sig_genes_cond,snps_genes) 

hdl_sig_genes_cond <- left_join(hdl_sig_genes_cond,
                                trans_ancestry_assign_ind[,c(1:2)],
                                by=c("rsID"="rsid_dbSNP150"))

hdl_sig_genes_cond <- hdl_sig_genes_cond %>%
  group_by(gene_id,Annot,chr,start,end,Mendelian,P,P_cond,rsID,pos,D) %>% summarise(trait = paste2(trait, collapse=","))
colnames(hdl_sig_genes_cond)[12] <- "glgc_trait"

hdl_sig_genes_cond <- hdl_sig_genes_cond[hdl_sig_genes_cond$P_cond < 0.05/length(unique(hdl_sig_genes_cond$gene_id)),]
hdl_sig_genes_cond <- arrange(hdl_sig_genes_cond,chr, start)

write_csv(hdl_sig_genes_cond,file = "/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/HDL.norm_sig_genes_cond.csv")

##############################

setwd("/rprojectnb2/pelosolab/yuxuan/lncRNA/")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/mendelian_genes.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/glgc/trans_ancestry_assign_ind.Rdata")

mendelian_genes$start_100k <- mendelian_genes$start - 100000
mendelian_genes$end_100k <- mendelian_genes$end + 100000

setDT(mendelian_genes)
setkey(mendelian_genes,chr,start_100k,end_100k)

path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results")
trait <- "TOTAL_ADJ.norm_cond"
annot <- c("gencode","fantom","lncrnakb","noncode")


gencode_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[1],".csv"))
fantom_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[2],".csv"))
#fantom_lncrna <- fantom_lncrna[-grep("ENSG", fantom_lncrna$gene_id),]
lncrnakb_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[3],".csv"))
noncode_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[4],".csv"))

lncrna_genes <- rbind(gencode_lncrna,fantom_lncrna,lncrnakb_lncrna,noncode_lncrna)
lncrna_genes <- lncrna_genes[,c("gene_id","STAAR-O")]

total_sig_genes <- read_csv("1_results/TOTAL_ADJ.norm_sig_genes.csv")
total_sig_genes <- total_sig_genes[,c(9,13,2,10:11,3,12)]
colnames(total_sig_genes) <- c("gene_id","Annot","chr","start","end","Mendelian","STAAR-O")

total_sig_genes <- left_join(total_sig_genes,lncrna_genes,by = c("gene_id"))

total_sig_genes_cond <- rbind(total_sig_genes[is.na(total_sig_genes$Mendelian),],
                              aggregate(Mendelian ~ gene_id + Annot + chr + start + end + `STAAR-O.x` + `STAAR-O.y`,
                                        data = total_sig_genes, function(x) paste2(x, collapse=";")))
total_sig_genes_cond <- arrange(total_sig_genes_cond,chr,start)
colnames(total_sig_genes_cond)[c(7,8)] <- c("P","P_cond")

genes_bed <- total_sig_genes_cond[,c("chr","start", "end", "gene_id")]
genes_bed$chr <- paste0("chr",genes_bed$chr)
write.table(genes_bed, "2_results/genes_all.bed", row.names = F,col.names = F,quote = F,sep = '\t')

load("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/lipid_known_loci_info_rsID.Rdata")
snps_bed <- known_loci_info[,c("CHR","POS","POS","rsID")]
snps_bed <- snps_bed[!duplicated(snps_bed), ]
snps_bed$CHR <- paste0("chr",snps_bed$CHR)
write.table(snps_bed , "2_results/snps_all.bed", row.names = F,col.names = F,quote = F,sep = '\t')

setwd("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results")
Sys.setenv(PATH="/usr/bin/")
system("sort -k1,1 -k2,2n genes_all.bed > genes_all_sort.bed")
system("sort -k1,1 -k2,2n snps_all.bed > snps_all_sort.bed")

Sys.setenv(PATH="/share/pkg.7/bedtools/2.30.0/install/bin/")
system("bedtools closest -a genes_all_sort.bed -b snps_all_sort.bed  -t all -d -k 1 > snps_genes.bed")
Sys.setenv(PATH="/usr/bin/")

snps_genes <- read.table(file = "snps_genes.bed")
snps_genes <- snps_genes[,c("V4","V8","V7","V9")]
colnames(snps_genes) <- c("gene_id","rsID","pos","D")
total_sig_genes_cond <- left_join(total_sig_genes_cond,snps_genes) 

total_sig_genes_cond <- left_join(total_sig_genes_cond,
                                  trans_ancestry_assign_ind[,c(1:2)],
                                  by=c("rsID"="rsid_dbSNP150"))

total_sig_genes_cond <- total_sig_genes_cond %>%
  group_by(gene_id,Annot,chr,start,end,Mendelian,P,P_cond,rsID,pos,D) %>% summarise(trait = paste2(trait, collapse=","))
colnames(total_sig_genes_cond)[12] <- "glgc_trait"

total_sig_genes_cond <- total_sig_genes_cond[total_sig_genes_cond$P_cond < 0.05/length(unique(total_sig_genes_cond$gene_id)),]
total_sig_genes_cond <- arrange(total_sig_genes_cond,chr, start)

write_csv(total_sig_genes_cond,file = "/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/TOTAL_ADJ.norm_sig_genes_cond.csv")

#################

setwd("/rprojectnb2/pelosolab/yuxuan/lncRNA/")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/mendelian_genes.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/glgc/trans_ancestry_assign_ind.Rdata")

mendelian_genes$start_100k <- mendelian_genes$start - 100000
mendelian_genes$end_100k <- mendelian_genes$end + 100000

setDT(mendelian_genes)
setkey(mendelian_genes,chr,start_100k,end_100k)

path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results")
trait <- "TG_LOG.norm_cond"
annot <- c("gencode","fantom","lncrnakb","noncode")


gencode_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[1],".csv"))
fantom_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[2],".csv"))
#fantom_lncrna <- fantom_lncrna[-grep("ENSG", fantom_lncrna$gene_id),]
lncrnakb_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[3],".csv"))
noncode_lncrna <- read_csv(file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",trait,"_",annot[4],".csv"))

lncrna_genes <- rbind(gencode_lncrna,fantom_lncrna,lncrnakb_lncrna,noncode_lncrna)
lncrna_genes <- lncrna_genes[,c("gene_id","STAAR-O")]

logtg_sig_genes <- read_csv("1_results/TG_LOG.norm_sig_genes.csv")
logtg_sig_genes <- logtg_sig_genes[,c(9,13,2,10:11,3,12)]
colnames(logtg_sig_genes) <- c("gene_id","Annot","chr","start","end","Mendelian","STAAR-O")

logtg_sig_genes <- left_join(logtg_sig_genes,lncrna_genes,by = c("gene_id"))

logtg_sig_genes_cond <- rbind(logtg_sig_genes[is.na(logtg_sig_genes$Mendelian),],
                              aggregate(Mendelian ~ gene_id + Annot + chr + start + end + `STAAR-O.x` + `STAAR-O.y`,
                                        data = logtg_sig_genes, function(x) paste2(x, collapse=";")))
logtg_sig_genes_cond <- arrange(logtg_sig_genes_cond,chr,start)
colnames(logtg_sig_genes_cond)[c(7,8)] <- c("P","P_cond")

genes_bed <- logtg_sig_genes_cond[,c("chr","start", "end", "gene_id")]
genes_bed$chr <- paste0("chr",genes_bed$chr)
write.table(genes_bed, "2_results/genes_all.bed", row.names = F,col.names = F,quote = F,sep = '\t')

load("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/lipid_known_loci_info_rsID.Rdata")
snps_bed <- known_loci_info[,c("CHR","POS","POS","rsID")]
snps_bed <- snps_bed[!duplicated(snps_bed), ]
snps_bed$CHR <- paste0("chr",snps_bed$CHR)
write.table(snps_bed , "2_results/snps_all.bed", row.names = F,col.names = F,quote = F,sep = '\t')

setwd("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results")
Sys.setenv(PATH="/usr/bin/")
system("sort -k1,1 -k2,2n genes_all.bed > genes_all_sort.bed")
system("sort -k1,1 -k2,2n snps_all.bed > snps_all_sort.bed")

Sys.setenv(PATH="/share/pkg.7/bedtools/2.30.0/install/bin/")
system("bedtools closest -a genes_all_sort.bed -b snps_all_sort.bed  -t all -d -k 1 > snps_genes.bed")
Sys.setenv(PATH="/usr/bin/")

snps_genes <- read.table(file = "snps_genes.bed")
snps_genes <- snps_genes[,c("V4","V8","V7","V9")]
colnames(snps_genes) <- c("gene_id","rsID","pos","D")
logtg_sig_genes_cond <- left_join(logtg_sig_genes_cond,snps_genes) 

logtg_sig_genes_cond <- left_join(logtg_sig_genes_cond,
                                  trans_ancestry_assign_ind[,c(1:2)],
                                  by=c("rsID"="rsid_dbSNP150"))

logtg_sig_genes_cond <- logtg_sig_genes_cond %>%
  group_by(gene_id,Annot,chr,start,end,Mendelian,P,P_cond,rsID,pos,D) %>% summarise(trait = paste2(trait, collapse=","))
colnames(logtg_sig_genes_cond)[12] <- "glgc_trait"

logtg_sig_genes_cond <- logtg_sig_genes_cond[logtg_sig_genes_cond$P_cond < 0.05/length(unique(logtg_sig_genes_cond$gene_id)),]
logtg_sig_genes_cond <- arrange(logtg_sig_genes_cond,chr, start)

write_csv(logtg_sig_genes_cond,file = "/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/TG_LOG.norm_sig_genes_cond.csv")
