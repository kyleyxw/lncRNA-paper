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
#pheno <- read_csv("/rprojectnb2/pelosolab/topmed/topmed_freeze.8.lipids.for_analysis.2020Feb11.csv")
lipid <- args$arg1 # "LDL_ADJ.norm", "HDL.norm", "TOTAL_ADJ.norm", "TG_LOG.norm", "nonHDL.norm","lTG_HDL.norm"
chr <- args$arg2 # chr1-22
annot <- args$arg3 #
type <- args$arg4

# lipid <- "LDL_ADJ.norm" # "LDL_ADJ.norm", "HDL.norm", "TOTAL_ADJ.norm", "TG_LOG.norm", "nonHDL.norm","lTG_HDL.norm"
# chr <- "chr4" # chr1-22
# annot <- "ensembl_gene" #
# type <- "gene"

print(paste(annot, lipid, chr,"is running"))


# ncRNAs by annot
lncrna <- get(load(paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/annot/","lncrna_",annot,"_",chr,".Rdata")))


# genofile
gdsfile <- paste0("/restricted/projectnb/pelosolab/topmed/frz8_agds/freeze.8.",chr,".pass_and_fail.gtonly.minDP0.gds")
genofile <- seqOpen(gdsfile)

# null model
obj_nullmodel <- readRDS(paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/null_model/",lipid, ".rds"))

phenotype.id <- as.character(obj_nullmodel$id_include)

# known_loci

known_loci <- get(load("/restricted/projectnb/pelosolab/yuxuan/lncRNA/2_results/lipid_known_loci_info.Rdata"))
known_loci$CHR <- paste0("chr",known_loci$CHR)

## get SNV id
filter <- seqGetData(genofile, "annotation/filter")
AVGDP <- seqGetData(genofile, "annotation/info/AVGDP")
SNVlist <- filter == "PASS" & AVGDP > 10 & isSNV(genofile)
variant.id <- seqGetData(genofile, "variant.id")

rm(AVGDP)
gc()

rm(filter)
gc()

## ncRNA SNVs
#GENCODE.Category <- seqGetData(genofile, "annotation/info/TOPMedAnnotation/GENCODE.Category")
#is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)

is.in <- SNVlist

variant.id <- variant.id[is.in]

#rm(GENCODE.Category)
#gc()

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

ncRNA_cond <- function(chr,gene_name,genofile,obj_nullmodel,
                       rare_maf_cutoff=0.01,rv_num_cutoff=2,
                       method_cond=c("optimal","naive"),
                       geno_missing_imputation=c("mean","minor")){
  
  seqResetFilter(genofile)
  
  ## evaluate choices
  method_cond <- match.arg(method_cond)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  position <- seqGetData(genofile, "position")
  REF <- as.character(seqGetData(genofile, "$ref"))
  ALT <- as.character(seqGetData(genofile, "$alt"))
  variant.id <- seqGetData(genofile, "variant.id")
 
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
  
  ## known loci
  known_loci_chr <- known_loci[known_loci[,1]==chr,]
  #known_loci_chr <- known_loci_chr[order(known_loci_chr[,2]),]
  
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
  
  ## Genotype Info
  REF_region <- as.character(seqGetData(genofile, "$ref"))
  ALT_region <- as.character(seqGetData(genofile, "$alt"))
  
  position_region <- as.numeric(seqGetData(genofile, "position"))
  sub_start_loc <- min(position_region)
  sub_end_loc <- max(position_region)
  
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
  results <- c()
  
  ### known variants needed to be adjusted
  known_loci_chr_region <- known_loci_chr[(known_loci_chr[,2]>=sub_start_loc-1E6)&(known_loci_chr[,2]<=sub_end_loc+1E6),]
  
  if(dim(known_loci_chr_region)[1]==0){
    
    results_temp <- c()
    results <- rbind(results,results_temp)
    
    }else
    {
      ## Genotype of Adjusted Variants
      rs_num_in <- c()
      for(i in 1:dim(known_loci_chr_region)[1])
      {
        rs_num_in <- c(rs_num_in,which((position==known_loci_chr_region[i,2])&(REF==known_loci_chr_region[i,3])&(ALT==known_loci_chr_region[i,4])))
      }
      
      variant.id.in <- variant.id[rs_num_in]
      seqSetFilter(genofile,variant.id=variant.id.in,sample.id=phenotype.id)
      
      ## genotype id
      id.genotype <- seqGetData(genofile,"sample.id")
      
      id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
      phenotype.id.merge <- data.frame(phenotype.id)
      phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
      id.genotype.match <- phenotype.id.merge$index
      
      Geno_adjusted <- seqGetData(genofile, "$dosage")
      Geno_adjusted <- Geno_adjusted[id.genotype.match,,drop=FALSE]
      
      ## impute missing
      if(!is.null(dim(Geno_adjusted)))
      {
        if(dim(Geno_adjusted)[2]>0)
        {
          if(geno_missing_imputation=="mean")
          {
            Geno_adjusted <- STAARpipeline::matrix_flip_mean(Geno_adjusted)$Geno
          }
          if(geno_missing_imputation=="minor")
          {
            Geno_adjusted <- STAARpipeline::matrix_flip_minor(Geno_adjusted)$Geno
          }
        }
      }
      
      if(class(Geno_adjusted)=="numeric")
      {
        Geno_adjusted <- matrix(Geno_adjusted,ncol=1)
      }
      
      AF <- apply(Geno_adjusted,2,mean)/2
      MAF <- AF*(AF<0.5) + (1-AF)*(AF>=0.5)
      
      Geno_adjusted <- Geno_adjusted[,MAF>0]
      if(class(Geno_adjusted)=="numeric")
      {
        Geno_adjusted <- matrix(Geno_adjusted,ncol=1)
      }
      
      
      seqResetFilter(genofile)
      
      ## Exclude RV in the region which needed to be adjusted
      id_exclude <- c()
      for(i in 1:length(rs_num_in))
      {
        id_exclude <- c(id_exclude,which((position_region==known_loci_chr_region[i,2])&(REF_region==known_loci_chr_region[i,3])&(ALT_region==known_loci_chr_region[i,4])))
      }
      
      if(length(id_exclude)>0)
      {
        Geno <- Geno[,-id_exclude]
        Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[-id_exclude,]
      }
      
      pvalues <- 0
      try(pvalues <- STAAR_cond(Geno,Geno_adjusted,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond=method_cond))
      
      
      if(class(pvalues)=="list")
      {
        results_temp <- rep(NA,4)
        results_temp[3] <- "lncRNA_cond"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        
        results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        
        results <- rbind(results,results_temp)
      }
      
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
  tmp <- ncRNA_cond(chr=chr,gene_name = genes[i],genofile = genofile,obj_nullmodel = obj_nullmodel,geno_missing_imputation="mean")
  results <- rbind(results, tmp)
  gc()
}

results <- apply(results,2,as.character)
results <- data.frame(results,check.names = F)


# Write the first data set in a new workbook
save(results, file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/",lipid,"_cond","/topmed_lncRNA_",annot,"_", lipid,"_",chr, "_cond.Rdata"))
showfile.gds(closeall=TRUE)

sessionInfo()

