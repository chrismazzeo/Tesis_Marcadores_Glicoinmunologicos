#clear console
cat("\014") 
#clear enviroment
rm(list=ls())

list.of.packages <- c("readr","BiocParallel","dplyr")
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
library("dplyr")

#***************************************************************
#todo #settings
alphaValue = 0.05

HTCountsSourcePath = "../resources/USCS/Gene Expression"
HTCountsDetinyPath = "../results/"



bd <- read_delim(paste0(HTCountsSourcePath,"/TcgaTargetGtex_gene_expected_count"), 
                                                 "\t", escape_double = FALSE, trim_ws = TRUE)


phenotype <- read_delim("../resources/USCS/USCS-TCGA-TARGET-GTEx-COAD-PHENOTYPES.csv",";", escape_double = FALSE, trim_ws = TRUE)
View(phenotype)

#hay 2 samples que no estan como columnas en bd
a = colnames(bd)
a = a[order(a)]
caca = subset(a, a %in% phenotype$samples)
COAD_BD =as.data.frame(bd[,caca])
row.names(COAD_BD) = bd$sample
write.csv(COAD_BD, file="COAD_TcgaTargetGtex_gene_expected_count.csv")

fixedPhenotype = subset(phenotype, phenotype$sample %in% colnames(COAD_BD))
fixedPhenotype = fixedPhenotype[order(fixedPhenotype$sample),]
write.csv(fixedPhenotype, file="FIXED_USCS_COAD_PHENOTYPE.csv")




