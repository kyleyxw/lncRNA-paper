rm(list = ls())

library(readr)
library(ggplot2)
LDL_ADJ_norm_lncrnakb <- read_csv("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/LDL_ADJ.norm_lncrnakb.csv")
LDL_ADJ_norm_gencode <- read_csv("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/LDL_ADJ.norm_gencode.csv")
LDL_ADJ_norm_fantom <- read_csv("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/LDL_ADJ.norm_fantom.csv")
LDL_ADJ_norm_noncode <- read_csv("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/LDL_ADJ.norm_noncode.csv")

HDL_norm_sig_genes <- read_csv("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/HDL.norm_sig_genes.csv")
LDL_ADJ_norm_sig_genes <- read_csv("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/LDL_ADJ.norm_sig_genes.csv")
TG_LOG_norm_sig_genes <- read_csv("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/TG_LOG.norm_sig_genes.csv")
TOTAL_ADJ_norm_sig_genes <- read_csv("/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/TOTAL_ADJ.norm_sig_genes.csv")
all_sig <- rbind(HDL_norm_sig_genes,LDL_ADJ_norm_sig_genes,TG_LOG_norm_sig_genes,TOTAL_ADJ_norm_sig_genes)
a <- table(all_sig$i.group)
b <- c(17741,11349,59633,78166)

df <- data.frame(a,b)
colnames(df) <- c("Annotation","No.Sig LncRNAs","No.Total LncRNAs Tested")
ggplot(df, aes(x=`No.Total LncRNAs Tested`,y=`No.Sig LncRNAs`, shape=Annotation, color = Annotation)) + 
  geom_point(size = 2) +
  theme_classic() +
  ylim(0,80) + 
  xlim(0,80000)

