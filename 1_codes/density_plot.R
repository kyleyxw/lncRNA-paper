# Exploratory analysis
rm(list = ls())
library(ggplot2)
library(GenomicRanges)
library(karyoploteR)

#########Density plot##############
myfiles <- list.files("annot","*gene_chr11.Rdata", full.names=T)
for (i in 1:length(myfiles)) {
  load(myfiles[i])
}

lncrna_ensembl_gene <- makeGRangesFromDataFrame(df = lncrna_ensembl_gene)
lncrna_gencode_gene <- makeGRangesFromDataFrame(df = lncrna_gencode_gene)
lncrna_fantom_gene <- makeGRangesFromDataFrame(df = lncrna_fantom_gene)
lncrna_lncrnakb_gene <- makeGRangesFromDataFrame(df = lncrna_lncrnakb_gene)
lncrna_noncode_gene <- makeGRangesFromDataFrame(df = lncrna_noncode_gene)

png(file = "chr11_density_plot.png", width = 1800, height = 2800, res = 300)
kp <- plotKaryotype(plot.type=2, chromosomes = "chr11")
kp <- kpPlotDensity(kp, r0=0, r1=0.2,data=lncrna_ensembl_gene)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.15)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.15)
kp <- kpPlotDensity(kp, r0=0.2, r1=0.35,data=lncrna_gencode_gene)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.2, r1=0.35)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0.2, r1=0.35)
kp <- kpPlotDensity(kp, r0=0.4, r1=0.55,data=lncrna_fantom_gene)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.4, r1=0.55)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0.4, r1=0.55)
kp <- kpPlotDensity(kp, r0=0.6, r1=0.75,data=lncrna_lncrnakb_gene)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.6, r1=0.75)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0.6, r1=0.75)
kp <- kpPlotDensity(kp, r0=0.8, r1=0.95,data=lncrna_noncode_gene)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.8, r1=0.95)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0.8, r1=0.95)
dev.off()

#########Density plot##############

#y <- c(781,798,1285,3017,4191,0,0,1838,1807,1805,9231,7282,1227,7342)
y <- c(781,798,1285,3017,4191)
#x <- c("Ensembl","GENCODE","FANTOM","LncRNAKB","NONCODE","BIGTranscriptome","MiTranscriptome",
#       "Ensembl","GENCODE","FANTOM","LncRNAKB","NONCODE","BIGTranscriptome","MiTranscriptome")
x <- c("Ensembl","GENCODE","FANTOM","LncRNAKB","NONCODE")
feature <- c(rep("Gene",5))
df <- data.frame(x,y,feature)
df$x <- factor(df$x, levels = c("Ensembl","GENCODE","FANTOM","LncRNAKB","NONCODE"))

ggplot(df, aes(x=x,y=y,fill=feature)) + 
  geom_bar(position="dodge", stat="identity") + xlab("Annotation") +ylab("Number")+
  theme_classic() + scale_fill_grey()
ggsave("chr11_bar_plot.png",width = 6, height = 4)


#####LDL#####################
source("/rprojectnb2/pelosolab/yuxuan/utilities/qq_plot.R")
myfiles <- list.files("chr11_v3","*gene_LDL_ADJ.norm_chr11.Rdata", full.names=T)
annot <- c("Ensembl","FANTOM","GENCODE","LncRNAKB","NONCODE")
for (i in 1:length(myfiles)) {
  load(myfiles[i])
  results[ , c(4:90)] <- apply(results[ ,c(4:90)], 2, function(x) as.numeric(x))
  my.pvalue.list<-list("All"=results$`STAAR-O`)
  qqplot <- paste("plot1_", i, sep = "")
  assign(qqplot, qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)),main=paste0(annot[i],"_lncRNA_gene_LDL_ADJ.norm_chr11")))
}

for (i in 1:length(myfiles)) {
  load(myfiles[i])
  results[ , c(4:90)] <- apply(results[ ,c(4:90)], 2, function(x) as.numeric(x))
  results$group <- cut(results$`#SNV`, breaks = c(0, 200, 1000, 5000))
  my.pvalue.list<-list("1-200"=results[results$group=="(0,200]",]$`STAAR-O`,
                       "201-1000"=results[results$group=="(200,1e+03]",]$`STAAR-O`,
                       "1001-5000"=results[results$group=="(1e+03,5e+03]",]$`STAAR-O`)
  qqplot <- paste("plot2_", i, sep = "")
  assign(qqplot, qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)),main=paste0(annot[i],"_lncRNA_gene_LDL_ADJ.norm_chr11")))
}

#################### Margaret########
path_to_ldl <- "/rprojectnb2/pelosolab/yuxuan/lncRNA/LDL_margaret/"

load(paste0(path_to_ldl,"ncRNA_62.Rdata"))
results <- apply(results_ncRNA,2,as.character)
ncrna_62 <- data.frame(results,check.names = F)

load(paste0(path_to_ldl,"ncRNA_63.Rdata"))
results <- apply(results_ncRNA,2,as.character)
ncrna_63 <- data.frame(results,check.names = F)

load(paste0(path_to_ldl,"ncRNA_64.Rdata"))
results <- apply(results_ncRNA,2,as.character)
ncrna_64 <- data.frame(results,check.names = F)

load(paste0(path_to_ldl,"ncRNA_65.Rdata"))
results <- apply(results_ncRNA,2,as.character)
ncrna_65 <- data.frame(results,check.names = F)

load(paste0(path_to_ldl,"ncRNA_66.Rdata"))
results <- apply(results_ncRNA,2,as.character)
ncrna_66 <- data.frame(results,check.names = F)

results <- rbind(ncrna_62,ncrna_63,ncrna_64,ncrna_65,ncrna_66)
results[ , c(4:90)] <- apply(results[ ,c(4:90)], 2, function(x) as.numeric(x))
results$group <- cut(results$`#SNV`, breaks = c(0, 200, 1000,5000))
my.pvalue.list<-list("All"=results$`STAAR-O`)
plot1_6 <- qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)),main="LDL_STAAR_Pipeline_ncRNA")
my.pvalue.list<-list("1-200"=results[results$group=="(0,200]",]$`STAAR-O`,
                     "201-1000"=results[results$group=="(200,1e+03]",]$`STAAR-O`,
                     "1001-5000"=results[results$group=="(1e+03,5e+03]",]$`STAAR-O`)
plot2_6 <- qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)),main="LDL_STAAR_Pipeline_ncRNA")


png("ldl_1.png", width = 1200, height = 800)
grid.arrange(plot1_1,plot1_2,plot1_3,plot1_4,plot1_5,plot1_6, ncol=3)
dev.off()

png("ldl_2.png", width = 1200, height = 800)
grid.arrange(plot2_1,plot2_2,plot2_3,plot2_4,plot2_5,plot2_6, ncol=3)
dev.off()
