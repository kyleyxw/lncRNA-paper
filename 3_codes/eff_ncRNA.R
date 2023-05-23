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


eff_ncRNA <- function(chr,gene_name,genofile,obj_nullmodel,
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
  
  if(class(Geno) != "matrix" && !(!is.null(attr(class(Geno), "package")) && attr(class(Geno), "package") == "Matrix")){
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

results  <- data.frame(matrix(ncol = 0, nrow = 63424))
for (i in 1:length(genes)){ 
  print(paste("No",i,"gene is running"))
  try(tmp <- eff_ncRNA(chr=chr,gene_name = genes[i],genofile = genofile,obj_nullmodel = obj_nullmodel,geno_missing_imputation="mean"))
  results <- cbind(results, tmp)
  gc()
}


save(results, file = paste0("/rprojectnb2/pelosolab/yuxuan/lncRNA/3_results/topmed_lncRNA_",annot,"_",chr, "_eff.Rdata"))
showfile.gds(closeall=TRUE)

sessionInfo()

