# samples
rm(list = ls())
library(readr)
library(GWASTools)
library(data.table)
library(Biobase)
library(SeqArray)
library(SeqVarTools)
library(STAAR)


lipid <- c( "LDL_ADJ.norm", "HDL.norm", "TOTAL_ADJ.norm", "TG_LOG.norm", "nonHDL.norm","lTG_HDL.norm") 
pheno <- read_csv("/topmed/topmed_freeze.8.lipids.for_analysis.2020Feb11.csv")
sample.id <- pheno$sample.id

#fit null model
## load GRM
load("pcrelate_kinshipMatrix_sparseDeg4.RData")
kmat <- km[sample.id,sample.id]

for (i in lipid){
  fomula <- as.formula(paste(i, " ~ ", paste(c("age","age2","sex","study",paste("PC", 1:11, sep="")), collapse= "+")))
  obj_nullmodel <- fit_null_glmmkin(fomula,data=pheno,family=gaussian(link = "identity"),id = "sample.id", kins = kmat)
  saveRDS(obj_nullmodel, file = paste0("lncRNA/null_model/",i, ".rds"))
}  

