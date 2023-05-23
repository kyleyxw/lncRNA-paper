rm(list = ls())
library(meta)
library(dplyr)
library(tidyverse)
library(qqman)
library(CMplot)
library(ggrepel)
library(ggman)
library(ggpubr)


load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_fantom_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_ensembl_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_gencode_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_lncrnakb_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_noncode_gene.Rdata")

lipid <- c("LDL_ADJ.norm", "HDL.norm", "TOTAL_ADJ.norm", "TG_LOG.norm")
for (i in lipid){
  path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/",i)
  setwd(path)
  
  # Ensembl
  ## Ensembl gene
  
  filelist <- list.files(path,pattern="topmed_lncRNA_ensembl_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  gwasResults <- inner_join(gwasResults, lncrna_ensembl_gene)
  gwasResults <- gwasResults[,c(1,2,91,90)]
  gwasResults$chr <- as.numeric(str_remove(gwasResults$chr, "chr"))
  gwasResults$`STAAR-O` <- as.numeric(as.character(gwasResults$`STAAR-O`))
  colnames(gwasResults) <- c("SNP","CHR","BP","P")
  snpsOfInterest <- gwasResults[gwasResults$P < (0.05 / nrow(gwasResults)),]
  p1 <- ggman(gwasResults, snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P",
              sigLine = -log10(0.05 / nrow(gwasResults)), lineColour = 'black',relative.positions = F, 
              title = paste0("Manhattan Plot ", i, " Ensembl LncRNAs")) + 
    #scale_color_manual(values = rep(c("grey30", "deepskyblue4"), 22 )) +
    scale_color_discrete()+
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  p1 <- ggmanLabel(p1, labelDfm = snpsOfInterest, snp = "SNP", label = "SNP",type = "text")
  
  ## Gencode gene
  
  filelist <- list.files(path,pattern="topmed_lncRNA_gencode_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  gwasResults <- inner_join(gwasResults, lncrna_gencode_gene)
  gwasResults <- gwasResults[,c(1,2,91,90)]
  gwasResults$chr <- as.numeric(str_remove(gwasResults$chr, "chr"))
  gwasResults$`STAAR-O` <- as.numeric(as.character(gwasResults$`STAAR-O`))
  colnames(gwasResults) <- c("SNP","CHR","BP","P")
  #subset only the SNPs 
  snpsOfInterest <- gwasResults[gwasResults$P < (0.05 / nrow(gwasResults)),]
  p2 <- ggman(gwasResults, snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P",
              sigLine = -log10(0.05 / nrow(gwasResults)), lineColour = 'black',relative.positions = F, 
              title = paste0("Manhattan Plot ", i, " Gencode LncRNAs")) + 
    #scale_color_manual(values = rep(c("grey30", "deepskyblue4"), 22 )) +
    scale_color_discrete()+
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  p2 <- ggmanLabel(p2, labelDfm = snpsOfInterest, snp = "SNP", label = "SNP",type = "text")
  
  ## Noncode gene
  
  filelist <- list.files(path,pattern="topmed_lncRNA_noncode_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  gwasResults <- inner_join(gwasResults, lncrna_noncode_gene)
  gwasResults <- gwasResults[,c(1,2,91,90)]
  gwasResults$chr <- as.numeric(str_remove(gwasResults$chr, "chr"))
  gwasResults$`STAAR-O` <- as.numeric(as.character(gwasResults$`STAAR-O`))
  
  colnames(gwasResults) <- c("SNP","CHR","BP","P")
  #subset only the SNPs 
  snpsOfInterest <- gwasResults[gwasResults$P < (0.05 / nrow(gwasResults)),]
  p3 <- ggman(gwasResults, snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P",
              sigLine = -log10(0.05 / nrow(gwasResults)), lineColour = 'black',relative.positions = F, 
              title = paste0("Manhattan Plot ", i, " Noncode LncRNAs")) + 
    #scale_color_manual(values = rep(c("grey30", "deepskyblue4"), 22 )) +
    scale_color_discrete()+
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  p3 <- ggmanLabel(p3, labelDfm = snpsOfInterest, snp = "SNP", label = "SNP",type = "text")
  
  
  ## lncrnakb gene
  
  filelist <- list.files(path,pattern="topmed_lncRNA_lncrnakb_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  gwasResults <- inner_join(gwasResults, lncrna_lncrnakb_gene)
  gwasResults <- gwasResults[,c(1,2,91,90)]
  gwasResults$chr <- as.numeric(str_remove(gwasResults$chr, "chr"))
  gwasResults$`STAAR-O` <- as.numeric(as.character(gwasResults$`STAAR-O`))
  
  colnames(gwasResults) <- c("SNP","CHR","BP","P")
  #subset only the SNPs 
  snpsOfInterest <- gwasResults[gwasResults$P < (0.05 / nrow(gwasResults)),]
  p4 <- ggman(gwasResults, snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P",
              sigLine = -log10(0.05 / nrow(gwasResults)), lineColour = 'black',relative.positions = F, 
              title = paste0("Manhattan Plot ", i, " LncrnaKB LncRNAs")) + 
    #scale_color_manual(values = rep(c("grey30", "deepskyblue4"), 22 )) +
    scale_color_discrete()+
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  p4 <- ggmanLabel(p4, labelDfm = snpsOfInterest, snp = "SNP", label = "SNP",type = "text")
  
  
  ## Fantom gene
  
  
  filelist <- list.files(path,pattern="topmed_lncRNA_fantom_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  gwasResults <- inner_join(gwasResults, lncrna_fantom_gene)
  gwasResults <- gwasResults[,c(1,2,91,90)]
  gwasResults$chr <- as.numeric(str_remove(gwasResults$chr, "chr"))
  gwasResults$`STAAR-O` <- as.numeric(as.character(gwasResults$`STAAR-O`))
  
  colnames(gwasResults) <- c("SNP","CHR","BP","P")
  #subset only the SNPs 
  snpsOfInterest <- gwasResults[gwasResults$P < (0.05 / nrow(gwasResults)),]
  p5 <- ggman(gwasResults, snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P",
              sigLine = -log10(0.05 / nrow(gwasResults)), lineColour = 'black',relative.positions = F, 
              title = paste0("Manhattan Plot ", i, " Fantom LncRNAs")) + 
    #scale_color_manual(values = rep(c("grey30", "deepskyblue4"), 22 )) +
    scale_color_discrete()+
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  p5 <- ggmanLabel(p5, labelDfm = snpsOfInterest, snp = "SNP", label = "SNP",type = "text")
  
  ggarrange(p1, p3, p2, p4, p5, ncol = 2, nrow = 3)
  ggsave(paste0("1_plots/",i,".jpg"),width = 18, height = 12)
  
}

###################
path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/LDL_margaret/")
filelist <- list.files(path,pattern=".Rdata",full.names = TRUE)
gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
gwasResults <- do.call(rbind, gwasResults)
unlist(gwasResults)
gwasResults <- as.data.frame(gwasResults)

gwasResults <- data.frame(apply(gwasResults,2,unlist))
gwasResults <- inner_join(gwasResults, lncrna_gencode_gene,by = c("Gene.name" = "gene_name"))
gwasResults <- gwasResults[,c(1,2,91,90)]
gwasResults$chr <- as.numeric(str_remove(gwasResults$Chr, "Chr"))
gwasResults <- gwasResults[,c(1,3,4,5)]
gwasResults$STAAR.O <- as.numeric(as.character(gwasResults$STAAR.O))

colnames(gwasResults) <- c("SNP","BP","P","CHR")
#subset only the SNPs 
snpsOfInterest <- gwasResults[gwasResults$P < (0.05 / nrow(gwasResults)),]
p6 <- ggman(gwasResults, snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P",
            sigLine = -log10(0.05 / nrow(gwasResults)), lineColour = 'black',relative.positions = F, 
            title = paste0("Manhattan Plot LDL Margaret Results")) + 
  #scale_color_manual(values = rep(c("grey30", "deepskyblue4"), 22 )) +
  scale_color_discrete()+
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

p6 <- ggmanLabel(p6, labelDfm = snpsOfInterest, snp = "SNP", label = "SNP",type = "text")
