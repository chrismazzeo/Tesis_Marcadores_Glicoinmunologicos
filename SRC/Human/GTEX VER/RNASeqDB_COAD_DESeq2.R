#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

#clear console
cat("\014") 
#clear enviroment
rm(list=ls())

list.of.packages <- c("readr","BiocParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocParallel")
#biocLite("apeglm")
#biocLite("ReportingTools")
#biocLite("org.Hs.eg.db")

library(readr)
library(BiocParallel)
library(apeglm)
library(psych)
library(pastecs)
library(ReportingTools)
library(AnnotationDbi)
#***************************************************************
# INPUT 
#***************************************************************
colon_rsem_count_gtex <- read_delim("/Volumes/Externo/Projects/R pipelines/PFI ROADMAP/Resources/RNAseqDB/colon-rsem-count-gtex.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
colon_rsem_count_gtex = as.data.frame (colon_rsem_count_gtex)
rownames(colon_rsem_count_gtex) = colon_rsem_count_gtex$Hugo_Symbol
colon_rsem_count_gtex = colon_rsem_count_gtex[,-1]
colon_rsem_count_gtex = colon_rsem_count_gtex[,-1]
colon_rsem_count_gtex[1:3,1:3]


coad_rsem_count_tcga_normal <- read_delim("../Resources/RNAseqDB/coad-rsem-count-tcga.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
coad_rsem_count_tcga_normal = as.data.frame (coad_rsem_count_tcga_normal)
rownames(coad_rsem_count_tcga_normal) = coad_rsem_count_tcga_normal$Hugo_Symbol
coad_rsem_count_tcga_normal = coad_rsem_count_tcga_normal[,-1]
coad_rsem_count_tcga_normal = coad_rsem_count_tcga_normal[,-1]
coad_rsem_count_tcga_normal[1:3,1:3]


coad_rsem_count_tcga_tumor <- read_delim("../Resources/RNAseqDB/coad-rsem-count-tcga-t.txt", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE)
coad_rsem_count_tcga_tumor = as.data.frame (coad_rsem_count_tcga_tumor)
rownames(coad_rsem_count_tcga_tumor) = coad_rsem_count_tcga_tumor$Hugo_Symbol
coad_rsem_count_tcga_tumor = coad_rsem_count_tcga_tumor[,-1]
coad_rsem_count_tcga_tumor = coad_rsem_count_tcga_tumor[,-1]
coad_rsem_count_tcga_tumor[1:3,1:3]

#filtrar por la lista, me quedare con los que haya 
RNASeqDB_COAD_PHENOTYPE <- read_csv("../Resources/rnaseqdb/RNASeqDB_COAD_TRANSVERSE_phenotypes.csv")

GTEX_BD_FILTERED =as.data.frame(colon_rsem_count_gtex[,RNASeqDB_COAD_PHENOTYPE$sample])
GTEX_BD_FILTERED[1:3,1:3]

bd = cbind(GTEX_BD_FILTERED,coad_rsem_count_tcga_normal,coad_rsem_count_tcga_tumor)
bd = round(bd,0)

GTEX_COAD <- vector(mode="character", length=dim(GTEX_BD_FILTERED)[2])
GTEX_COAD[] = "GTEX_COAD"

TCGA_NORMAL <- vector(mode="character", length=dim(coad_rsem_count_tcga_normal)[2])
TCGA_NORMAL[] = "TCGA_COAD_NORMAL"

TCGA_TUMOR <- vector(mode="character", length=dim(coad_rsem_count_tcga_tumor)[2])
TCGA_TUMOR[] = "TCGA_COAD_TUMOR"
group = data.frame (condition = c(GTEX_COAD,TCGA_NORMAL,TCGA_TUMOR))


#*******************************************************

sampleTable <- data.frame(sampleName = colnames(bd),
                          fileName =colnames(bd),
                          condition = group$condition)

dds = DESeqDataSetFromMatrix(countData = bd, 
                             colData =sampleTable,
                             design= ~condition)

#***************************************************************
#DF analisys
#***************************************************************

register(MulticoreParam(4))
dds <- DESeq(dds, parallel = TRUE)
alphaValue = 0.05
#***************************************************************
#Comparo de a 2 grupos
#***************************************************************
resTCGA_NORMAL_VS_GTEX <- results(dds, parallel = TRUE, pAdjustMethod = "BH", alpha = alphaValue,contrast=c("condition","TCGA_COAD_NORMAL","GTEX_COAD"))
resTCGA_TUMOR_VS_TCGA_NORMAL <- results(dds, parallel = TRUE, pAdjustMethod = "BH", alpha = alphaValue,contrast=c("condition","TCGA_COAD_TUMOR","TCGA_COAD_NORMAL"))
resTCGA_TUMORL_VS_GTEX <- results(dds, parallel = TRUE, pAdjustMethod = "BH", alpha = alphaValue,contrast=c("condition","TCGA_COAD_TUMOR","GTEX_COAD"))

resTCGA_NORMAL_VS_GTEX
resTCGA_TUMOR_VS_TCGA_NORMAL
resTCGA_TUMORL_VS_GTEX

#***************************************************************
#summary results
#***************************************************************
mcols(resTCGA_NORMAL_VS_GTEX)$description
summary(resTCGA_NORMAL_VS_GTEX)
sum(resTCGA_NORMAL_VS_GTEX$padj < alphaValue, na.rm=TRUE)
plotMA(resTCGA_NORMAL_VS_GTEX, ylim=c(-2,2))

mcols(resTCGA_TUMOR_VS_TCGA_NORMAL)$description
summary(resTCGA_TUMOR_VS_TCGA_NORMAL)
sum(resTCGA_TUMOR_VS_TCGA_NORMAL$padj < alphaValue, na.rm=TRUE)
plotMA(resTCGA_TUMOR_VS_TCGA_NORMAL, ylim=c(-2,2))

mcols(resTCGA_TUMORL_VS_GTEX)$description
summary(resTCGA_TUMORL_VS_GTEX)
sum(resTCGA_TUMORL_VS_GTEX$padj < alphaValue, na.rm=TRUE)
plotMA(resTCGA_TUMORL_VS_GTEX, ylim=c(-2,2))

#***************************************************************
#Guardamos
#***************************************************************
#guardamos count matrix
write.csv(assay(dds),file="../Results/DE/rnaseqdb/RNASeqDB_COUNT_MATRIX.csv")
#guardamos DE
write.csv(resTCGA_NORMAL_VS_GTEX,file=paste0("../Results/DE/rnaseqdb/","alpha",alphaValue,"/","TCGA_NORMAL_VS_GTEXT_LFC_alpha",alphaValue,".csv"))
write.csv(resTCGA_TUMOR_VS_TCGA_NORMAL,file=paste0("../Results/DE/rnaseqdb/","alpha",alphaValue,"/","TCGA_TUMOR_VS_TCGA_NORMAL_LFC_alpha",alphaValue,".csv"))
write.csv(resTCGA_TUMORL_VS_GTEX,file=paste0("../Results/DE/rnaseqdb/","alpha",alphaValue,"/","TCGA_TUMORL_VS_GTEX_LFC_alpha",alphaValue,".csv"))

