###########################################################
# Extract CHR, POS, REF and ALT from GLGC #rs 
# Yuxuan Wang
# Initiate date: 03/28/2022
# Current date: 03/10/2022
###########################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

###########################################################
#           User Input
###########################################################
## aGDS directory


## Input GWASCatalog variants (#rs)
load("/restricted/projectnb/pelosolab/yuxuan/glgc/trans_ancestry_assign_ind.Rdata")
known_loci <- unique(na.omit(trans_ancestry_assign_ind$rsid_dbSNP150))
## rs channel name in aGDS
rs_channel <- "annotation/info/TOPMedAnnotation/dbSNP_rs_num"
## output path
output_path <- "/restricted/projectnb/pelosolab/yuxuan/lncRNA/2_results/"
## output file name
output_file_name <- "lipid_known_loci_info_rsID"

###########################################################
#           Main Function 
###########################################################
known_loci_info <- c()

for(chr in 1:22)
{
  print(chr)
  
  print(paste("Chromosome:",chr))
  gds.path <- paste0("/restricted/projectnb/pelosolab/topmed/frz8_agds/freeze.8.chr",chr,".pass_and_fail.gtonly.minDP0.gds")
  genofile <- seqOpen(gds.path)
  
  variant.id <- seqGetData(genofile, "variant.id")
  rs_num <- seqGetData(genofile,rs_channel)
  
  rs_num_in <- rs_num%in%known_loci
  
  if(sum(rs_num_in) > 0)
  {
    variant.id.in <- variant.id[rs_num_in]
    
    rm(rs_num)
    gc()
    
    rm(variant.id)
    gc()
    
    seqSetFilter(genofile,variant.id=variant.id.in)
    
    ### Basic Info of Significant Loci
    rsID <- seqGetData(genofile,rs_channel)
    position <- as.numeric(seqGetData(genofile, "position"))
    REF <- as.character(seqGetData(genofile, "$ref"))
    ALT <- as.character(seqGetData(genofile, "$alt"))
    
    known_loci_info_chr <- data.frame(rsID =rsID, CHR=rep(chr,length(position)),POS=position,REF=REF,ALT=ALT)
    known_loci_info <- rbind(known_loci_info,known_loci_info_chr)	
  }
  seqClose(genofile)
}

# Output Info of GWASCatalog SNVs
save(known_loci_info,file=paste0(output_path,output_file_name,".Rdata"))
