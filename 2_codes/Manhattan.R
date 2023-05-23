rm(list = ls())
library(ghibli)
library(meta)
library(dplyr)
library(tidyverse)
library(qqman)
library(CMplot)
library(ggrepel)
library(ggman)
library(ggpubr)



lipid <- c("LDL_ADJ.norm", "HDL.norm", "TOTAL_ADJ.norm", "TG_LOG.norm")
lipid_cond <- c("LDL_ADJ.norm_cond", "HDL.norm_cond", "TOTAL_ADJ.norm_cond", "TG_LOG.norm_cond")

annot <- c("gencode","fantom","lncrnakb","noncode")

path <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/1_results/"
path_cond <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/"

load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_fantom_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_ensembl_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_gencode_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_lncrnakb_gene.Rdata")
load("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/lncrna_noncode_gene.Rdata")


for (i in lipid){
  for (j in lipid_cond){
    
    df <- read_csv(file = paste0(path,i,"_",annot[1],".csv"))
    gwasResults <- inner_join(df, lncrna_gencode_gene)
    gwasResults_gencode <- gwasResults[,c(1,2,91,90)] 
    
    df <- read_csv(file = paste0(path,i,"_",annot[2],".csv"))
    gwasResults <- inner_join(df, lncrna_fantom_gene)
    gwasResults_fantom <- gwasResults[,c(1,2,91,90)]
    
    df <- read_csv(file = paste0(path,i,"_",annot[3],".csv"))
    gwasResults <- inner_join(df, lncrna_lncrnakb_gene)
    gwasResults_lncrnakb <- gwasResults[,c(1,2,91,90)]
    
    df <- read_csv(file = paste0(path,i,"_",annot[4],".csv"))
    gwasResults <- inner_join(df, lncrna_noncode_gene)
    gwasResults_noncode <- gwasResults[,c(1,2,91,90)]
    rm(df)
    
    gwasResults <- rbind(gwasResults_gencode,gwasResults_fantom,gwasResults_lncrnakb,gwasResults_noncode)
    gwasResults$chr <- as.numeric(str_remove(gwasResults$chr, "chr"))
    gwasResults$`STAAR-O` <- as.numeric(as.character(gwasResults$`STAAR-O`))
    colnames(gwasResults) <- c("SNP","CHR","BP","P")
    snpsOfInterest <- gwasResults[gwasResults$P <  0.05/113587,]
    
    
    df <- read_csv(file = paste0(path_cond,j,"_",annot[1],".csv"))
    gwasResults_cond <- inner_join(df, lncrna_gencode_gene)
    gwasResults_cond_gencode <- gwasResults_cond[,c(1,2,91,90)] 
    
    df <- read_csv(file = paste0(path_cond,j,"_",annot[2],".csv"))
    gwasResults_cond <- inner_join(df, lncrna_fantom_gene)
    gwasResults_cond_fantom <- gwasResults_cond[,c(1,2,91,90)]
    
    df <- read_csv(file = paste0(path_cond,j,"_",annot[3],".csv"))
    gwasResults_cond <- inner_join(df, lncrna_lncrnakb_gene)
    gwasResults_cond_lncrnakb <- gwasResults_cond[,c(1,2,91,90)]
    
    df <- read_csv(file = paste0(path_cond,j,"_",annot[4],".csv"))
    gwasResults_cond <- inner_join(df, lncrna_noncode_gene)
    gwasResults_cond_noncode <- gwasResults_cond[,c(1,2,91,90)]
    rm(df)
    
    gwasResults_cond <- rbind(gwasResults_cond_gencode,gwasResults_cond_fantom,gwasResults_cond_lncrnakb,gwasResults_cond_noncode)
    gwasResults_cond$chr <- as.numeric(str_remove(gwasResults_cond$chr, "chr"))
    gwasResults_cond$`STAAR-O` <- as.numeric(as.character(gwasResults_cond$`STAAR-O`))
    colnames(gwasResults_cond) <- c("SNP","CHR","BP","P")
    
    snpsOfInterest_cond <- gwasResults_cond[gwasResults_cond$SNP %in% snpsOfInterest$SNP,]
    snpsOfInterest_cond$SNP <- paste0(snpsOfInterest_cond$SNP," ")
    snpsOfInterest_cond$BP <- snpsOfInterest_cond$BP + 100
    gwasResults <- data.frame(rbind(gwasResults,snpsOfInterest_cond))
    #gwasResults <- arrange(gwasResults,CHR,BP)
    
    
    snpsOfInterest$group <- "P"
    snpsOfInterest_cond$group <- "P_cond"
    
    
    snps <- rbind(snpsOfInterest,snpsOfInterest_cond)
    snps <- arrange(snps,CHR,BP)
    }
    
    set.seed(234)
    p1 <- ggman(gwasResults, snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P",
                sigLine = -log10(0.05 / nrow(gwasResults)), 
                lineColour = 'black',relative.positions = F,pointSize=0.5,
                title = paste0("Association of ",i," and LncRNAs")) +
      theme_bw() + 
      theme(
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
    
    highlightDfm <- p1[[1]]
    highlightDfm <- highlightDfm[highlightDfm$SNP %in% snps$SNP,1:4]
    highlightDfm <- left_join(highlightDfm,snps[,1:5])
    
    
    p1 <- ggmanHighlightGroup(p1, highlightDfm = highlightDfm, snp = "SNP", group = "group", 
                              legend.title = "Analysis Type",size = 2)
    
    p1 <- ggmanLabel(p1, labelDfm = highlightDfm, snp = "SNP", label = "SNP",type = "label",max.overlaps=7,
                     nudge_x = .15,
                     box.padding = 0.5,
                     nudge_y = 1,
                     segment.curvature = -0.1,
                     segment.ncp = 3,
                     segment.angle = 20) + 
      scale_color_manual(values = rep(c("grey30", "grey60"), 11))
    
    ggsave(plot = p1,filename = paste0("lncRNAs_STAAR_O_",i,".png"),width = 20,height = 12,device = "png")
    
}




#######################################


