library(data.table)
library(RSpectra)

path <- "lncRNA/"
annot <- c("gencode","fantom","noncode","lncrnakb")

for (chr in 22:1){
  score <- data.frame(seq(1,63424))
  
  for (i in annot){
    print(paste0(path,"topmed_lncRNA_",i,"_gene_chr",chr,"_eff.Rdata"))
    genes <-  get(load(paste0(path,"topmed_lncRNA_",i,"_gene_chr",chr,"_eff.Rdata")))
    tmp <- read.csv(paste0("lncRNA/LDL_ADJ.norm_",i,".csv"))
    tmp <- tmp[tmp$chr == paste0("chr",chr),]
    genes <- genes[,tmp$gene_id]
    score<-cbind(score,genes)
  }
  
  score<-score[,-1]
  corscore <-  WGCNA::cor(score,nThreads = 24, verbose = T) 
  fwrite(corscore,
         file = paste0("lncRNA/topmed_lncRNA_",chr, "_gene_cor.csv"),
         col.names = T,row.names = T)
  
  rm(genes)
  gc()
  
  # PCA
  eigenscore<-eigs_sym(A=corscore,k=ncol(score),which="LM")
  rm(corscore)
  gc()
  
  eigenscore<-eigenscore$values
  
  
  j <- 1  
  
  while ((sum(eigenscore[1:j])/sum(eigenscore[1:length(eigenscore)])) < 0.995) {
    j <- j + 1
  }
  
  print(paste0("chr", chr," No.Eff Genes=", j))
  
}























