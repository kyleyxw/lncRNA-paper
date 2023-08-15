rm(list = ls())
library(data.table)
library(readr)
library(SeqArray)
library(SeqVarTools)
library(STAAR)

args = commandArgs(TRUE)


## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL
rm(argsL)

# phenofile

lipid <- args$arg1 # "LDL_ADJ.norm", "HDL.norm", "TOTAL_ADJ.norm", "TG_LOG.norm", "nonHDL.norm","lTG_HDL.norm"
chr <- args$arg2 # chr1-22
annot <- args$arg3 #
type <- args$arg4

print(paste(annot, lipid, chr,"is running"))


# lncRNAs by annot
lncrna <- get(load(paste0("lncRNA/annot/","lncrna_",annot,"_",chr,".Rdata")))

# genofile
gdsfile <- paste0("topmed/frz8_agds/freeze.8.",chr,".pass_and_fail.gtonly.minDP0.gds")
genofile <- seqOpen(gdsfile)

# null model
obj_nullmodel <- readRDS(paste0("lncRNA/null_model/",lipid, ".rds"))

phenotype.id <- as.character(obj_nullmodel$id_include)

## get SNV id
filter <- seqGetData(genofile, "annotation/filter")
AVGDP <- seqGetData(genofile, "annotation/info/AVGDP")
SNVlist <- filter == "PASS" & AVGDP > 10 & isSNV(genofile)
variant.id <- seqGetData(genofile, "variant.id")

rm(AVGDP)
gc()

rm(filter)
gc()

is.in <- SNVlist

variant.id <- variant.id[is.in]

position <- seqGetData(genofile, "position")
position <- position[is.in]

snps_info <- data.table(start=position,end=position, variant.id=variant.id)

if (type == "gene") {
  genes_info <- data.table(start=lncrna$start,end=lncrna$end,gene_id=lncrna$gene_id)
  setkey(genes_info, start, end)
  overlap <- foverlaps(snps_info, genes_info, type="within",nomatch = 0)
  freq <- as.data.frame(table(overlap$gene_id))
  freq <- freq[freq$Freq < 5000,]
  genes <- as.character(freq$Var1)
} else {
  genes_info <- data.table(start=lncrna$start,end=lncrna$end,transcript_id=lncrna$transcript_id)
  setkey(genes_info, start, end)
  overlap <- foverlaps(snps_info, genes_info, type="within",nomatch = 0)
  freq <- as.data.frame(table(overlap$transcript_id))
  freq <- freq[freq$Freq < 5000,]
  genes <- as.character(freq$Var1)
}


rm(is.in)
gc()

lncRNA <- function(chr,gene_name,genofile,obj_nullmodel,
                  rare_maf_cutoff=0.01,rv_num_cutoff=2,geno_missing_imputation=c("mean","minor")){
  
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  ### Gene
  
  variant.is.in <- overlap[gene_id==gene_name]$variant.id

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
  
  ## Annotation
  CADD.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/CADD.FULL/PHRED")
  CADD.PHRED[is.na(CADD.PHRED)] <- 0
  LINSIGHT.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/LINSIGHT.PHRED.rounded")
  FATHMM.XF.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/FATHMM.XF.PHRED.rounded")
  APC1.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.EpigeneticActive")
  APC2.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.EpigeneticRepressed")
  APC3.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.EpigeneticTranscription")
  APC4.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.Conservation.v2")
  APC5.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.LocalDiversity.v2")
  APC6.PHRED <- -10*log10(1-10^(-APC5.PHRED/10))
  APC7.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.Mappability")
  APC8.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.TF")
  APC9.PHRED <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.Protein")
  
  Anno.Int.PHRED.sub <- data.frame(CADD.PHRED,LINSIGHT.PHRED,FATHMM.XF.PHRED,
                                   APC1.PHRED,APC2.PHRED,APC3.PHRED,APC4.PHRED,APC5.PHRED,
                                   APC6.PHRED,APC7.PHRED,APC8.PHRED,APC9.PHRED)
  
  pvalues <- 0
  try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff))
  
  results <- c()
  if(class(pvalues)=="list")
  {
    results_temp <- rep(NA,4)
    results_temp[3] <- "lncRNA"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    
    results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                      pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                      pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    
    results <- rbind(results,results_temp)
  }
  
  if(!is.null(results))
  {
    colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
    colnames(results)[1:4] <- c("gene_id","chr","category","number_variant")
    colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
  }
  
  seqResetFilter(genofile)
  return(results)
}

results  <- c()
for (i in 1:length(genes)){ 
  print(paste("No",i,"gene is running"))
  tmp <- lncRNA(chr=chr,gene_name = genes[i],genofile = genofile,obj_nullmodel = obj_nullmodel,geno_missing_imputation="mean")
  results <- rbind(results, tmp)
  gc()
}

results <- apply(results,2,as.character)
results <- data.frame(results,check.names = F)


# Write the first data set in a new workbook
save(results, file = paste0("lncRNA/",lipid,"/topmed_lncRNA_",annot,"_", lipid,"_",chr, ".Rdata"))
showfile.gds(closeall=TRUE)

sessionInfo()

