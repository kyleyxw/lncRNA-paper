## Non Synonymous Scores for lipid mendelian genes
## Material for null model with more covariates

rm(list = ls())
library(data.table)
library(readr)
library(dplyr)
library(SeqArray)
library(SeqVarTools)
library(STAAR)

score <- function(chr,gene_name,genofile,
                      rare_maf_cutoff=0.01,rv_num_cutoff=2,geno_missing_imputation=c("mean","minor")){
  
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  ### Gene
  
  variant.is.in <- overlap[gene_name==gene_name]$variant.id
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,]
  
  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- STAARpipeline::matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- STAARpipeline::matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  if(class(Geno)[1] != "matrix" && !(!is.null(attr(class(Geno), "package")) && attr(class(Geno), "package") == "Matrix")){
    stop("Geno is not a matrix!")
  }
  
  if(dim(Geno)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }
  
  if(!is.null(attr(class(Geno), "package")) && attr(class(Geno), "package") == "Matrix"){
    Geno <- as.matrix(Geno)
  }
  Geno <- matrix_flip(Geno)
  MAF <- Geno$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- Geno$Geno[,RV_label]
  
  paste0("number of variants is ", sum(RV_label))
  rm(Geno)
  gc()
  
  if(sum(RV_label) >= rv_num_cutoff){
    
    results <- data.frame(rowSums(Geno_rare))
    colnames(results) <- gene_name
    rownames(results) <- phenotype.id.merge$phenotype.id
    
    rm(Geno_rare)
    gc()
    
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }
  
  seqResetFilter(genofile)
  return(results)
}

closest_gene <- read.table(file = "closest_genes_to_lncRNA.bed")
closest_gene <- closest_gene %>% select(V4,V1,V2,V3,V9,V10,V7,V8,V11) %>% unique()
colnames(closest_gene) <- c("lncRNA","chr","lncRNA_start","lncRNA_end","gene_id","gene_name","gene_start","gene_end","distance")
closest_gene <- closest_gene[closest_gene$distance <= 250000,]

lipid_genes <- c("PCSK9","APOB","LPA","LPL","APOA5","APOC3","APOA1","PCSK7",
                 "CETP","PLA2G15","LDLR","ANGPTL8","TM6SF2","APOE","APOC1")
closest_gene$lipid_gene <- ifelse(closest_gene$gene_name %in% lipid_genes,1,0)
closest_gene <- closest_gene %>% group_by(lncRNA) %>% mutate(closest_rank = row_number())
closest_gene <- closest_gene[closest_gene$lipid_gene == 1 | closest_gene$closest_rank == 1, ]

gene_list <- unique(closest_gene[,c(2,5:8)])
colnames(gene_list)[c(4,5)] <- c("start","end")

results  <- data.frame(matrix(ncol = 0, nrow = 66329))

for (j in 1:nrow(gene_list)){
  
  chr <- as.character(gene_list$chr[j])
  gene_name <- as.character(gene_list$gene_name[j])
  
  # genofile
  gdsfile <- paste0("topmed/frz8_agds/freeze.8.",chr,".pass_and_fail.gtonly.minDP0.gds")
  genofile <- seqOpen(gdsfile)
  
  # null model
  pheno <- read_csv("topmed/topmed_freeze.8.lipids.for_analysis.2020Feb11.csv")
  
  phenotype.id <- as.character(pheno$sample.id)
  
  
  ## get SNV id
  filter <- seqGetData(genofile, "annotation/filter")
  AVGDP <- seqGetData(genofile, "annotation/info/AVGDP")
  SNVlist <- filter == "PASS" & AVGDP > 10 & isSNV(genofile)
  variant.id <- seqGetData(genofile, "variant.id")
  
  rm(AVGDP)
  gc()
  
  rm(filter)
  gc()
  
  ## nonsynonymous SNV
  GENCODE.EXONIC.Category  <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/GENCODE.EXONIC.Category")
  is.in <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(SNVlist)
  
  #is.in <- SNVlist
  
  variant.id <- variant.id[is.in]
  
  rm(GENCODE.EXONIC.Category)
  gc()
  
  position <- seqGetData(genofile, "position")
  position <- position[is.in]
  
  snps_info <- data.table(start=position,end=position, variant.id=variant.id)
  
  gene_coding <- gene_list[gene_list$gene_name == gene_name,]
  genes_info <- data.table(start=gene_coding$start,end=gene_coding$end,gene_name=gene_coding$gene_name)
  setkey(genes_info, start, end)
  overlap <- foverlaps(snps_info, genes_info, type="within",nomatch = 0)
  genes <- as.character(gene_coding$gene_name)
  
  
  rm(is.in)
  gc()
  
  
  try(tmp <- score(chr=chr,gene_name = genes,genofile = genofile, geno_missing_imputation="mean"))
  results <- cbind(results, tmp)
  gc()
  
  
  showfile.gds(closeall=TRUE)
  gc()
  
}

identical(pheno$sample.id,rownames(results))

pheno <- cbind(pheno,results)
write.csv(pheno,"lncRNA/null_model_non_synonymous/topmed_freeze.8.lipids.for_analysis.2023May22.csv",
          row.names = F,quote = F)
