library(readr)
library(tidyverse)

HDL_norm_sig_genes_cond <- read_csv("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/HDL.norm_sig_genes_cond.csv")
cor_ch11 <- data.table::fread("/rprojectnb2/pelosolab/yuxuan/lncRNA/3_results/topmed_lncRNA_11_gene_cor.csv")
library(tidyverse)
cor_ch11 <- cor_ch11 %>%
  remove_rownames() %>%
  column_to_rownames(var = 'V1')
cor_ch11 <- as.matrix(cor_ch11)

cor_ch11 <- cor_ch11[HDL_norm_sig_genes_cond[2:11,]$gene_id,HDL_norm_sig_genes_cond[2:11,]$gene_id]
