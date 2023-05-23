rm(list = ls())
library(readr)
lipid <- c("LDL_ADJ.norm_cond", "HDL.norm_cond", "TOTAL_ADJ.norm_cond", "TG_LOG.norm_cond")
#lipid <- c("LDL_ADJ.norm_cond")

for (i in lipid){
  path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/",i)
  
  filelist <- list.files(path,pattern="topmed_lncRNA_ensembl_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  
  write_csv(gwasResults,paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",i,"_ensembl.csv"))
  }

for (i in lipid){
  path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/",i)
  
  filelist <- list.files(path,pattern="topmed_lncRNA_gencode_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  
  write_csv(gwasResults,paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",i,"_gencode.csv"))
}

for (i in lipid){
  path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/",i)
  
  filelist <- list.files(path,pattern="topmed_lncRNA_fantom_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  
  write_csv(gwasResults,paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",i,"_fantom.csv"))
}

for (i in lipid){
  path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/",i)
  
  filelist <- list.files(path,pattern="topmed_lncRNA_lncrnakb_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  
  write_csv(gwasResults,paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",i,"_lncrnakb.csv"))
}

for (i in lipid){
  path <- paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/",i)
  
  filelist <- list.files(path,pattern="topmed_lncRNA_noncode_gene_",full.names = TRUE)
  gwasResults <- sapply(filelist, function(x) mget(load(x)), simplify = TRUE) 
  gwasResults <- do.call(rbind, gwasResults)
  
  write_csv(gwasResults,paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/2_results/",i,"_noncode.csv"))
}