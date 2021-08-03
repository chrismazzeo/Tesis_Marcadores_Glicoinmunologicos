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
library(tximport)
#***************************************************************
#todo #settings
alphaValue = 0.05


HTCountsDetinyPath = "../results/USCS"


print("****** Frecuencia *****")
freq<-table(samples$`Sample Type`)
prop<-round(prop.table(freq),2)
print( cbind(freq,prop))
cat(" \n")


COAD_TcgaTargetGtex_gene_expected_count <- read_csv("../Resources/USCS/Gene Expression/COAD_TcgaTargetGtex_gene_expected_count.csv")
COAD_TcgaTargetGtex_gene_expected_count = as.data.frame(COAD_TcgaTargetGtex_gene_expected_count)
row.names(COAD_TcgaTargetGtex_gene_expected_count) = COAD_TcgaTargetGtex_gene_expected_count$X1
COAD_TcgaTargetGtex_gene_expected_count = COAD_TcgaTargetGtex_gene_expected_count[,-1]
COAD_TcgaTargetGtex_gene_expected_count= COAD_TcgaTargetGtex_gene_expected_count[order(row.names(COAD_TcgaTargetGtex_gene_expected_count)),]




USCS_TCGA_TARGET_GTEx_COAD_PHENOTYPES <- read_csv("../Resources/USCS/FIXED_USCS_COAD_PHENOTYPE.csv")
USCS_TCGA_TARGET_GTEx_COAD_PHENOTYPES= as.data.frame(USCS_TCGA_TARGET_GTEx_COAD_PHENOTYPES)
USCS_TCGA_TARGET_GTEx_COAD_PHENOTYPES = USCS_TCGA_TARGET_GTEx_COAD_PHENOTYPES[,-1]
row.names(USCS_TCGA_TARGET_GTEx_COAD_PHENOTYPES) = USCS_TCGA_TARGET_GTEx_COAD_PHENOTYPES$sample

#***************************************************************
#Differential Expression con DESeq2
#Workflow http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#using-parallelization
#***************************************************************
# INPUT 
#La matriz que usemos para los reads counts tiene que ser un-normalized counts o estimated counts (si es signgle-end RNA-seq)
#O sea valores Enteros (nada de log)El valor de la fija j y columna i es el count, que es la cantidad de reads que se ele asigno a gen para esa muestra

#Count Matrix -->DESeqDataSetFromMatrix 
#Transcript quantification --> DESeqDataSetFromTximport (es el output de Salmon, Sailfish, kallisto o RSEM, tximport  uses estimated gene counts from the transcript abundance quantifiers, but not normalized counts.)
# HTsec-counts  --->DESeqDataSetFromHTSeq
# RangesSummarizedExperiment -->DESeqDataSet
#***************************************************************
library("DESeq2")
#genera una matriz donde las filas son los geneID, las columnas el nombre de las mustras, y adentro tiene el valor de los count 
countMatrix = as.matrix(COAD_TcgaTargetGtex_gene_expected_count)
dds = DESeqDataSetFromMatrix(countData = countMatrix, 
                             colData =USCS_TCGA_TARGET_GTEx_COAD_PHENOTYPES,
                             design= ~ `_sample_type`)


dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = HTCountsTempFolder,
                                  design= ~ condition)

#seteamos el grupo control
dds$condition <- relevel(dds$condition, ref = "Solid Tissue Normal")

#borro la carpeta temporal
unlink(HTCountsTempFolder, recursive = TRUE)


#Multifactor design
#si tenemos varios factores ejmplo single end y pairend end + treated y untreated
#donde type es un factor y condition es otro factor
#design(dds) <- formula(~ type + condition)
#***************************************************************
#DF analisys
#***************************************************************
#utilizo paralalelismo
register(MulticoreParam(4))

#hago el analisis de expresion diferencial
#This very simple function call does all the hard work. Briefly, this function performs three things:
#Compute a scaling factor for each sample to account for differences in read depth and complexity between samples
#Estimate the variance among samples
#Test for differences in expression among groups (time points, in our case)

dds <- DESeq(dds, parallel = TRUE)

#***************************************************************
#Comparo de a 2 grupos
#***************************************************************
#compare 2 groups
#seteamos:
#método de correcion por comparacion multiple
# alpha
# que grupos comparamos


res <- results(dds, parallel = TRUE, pAdjustMethod = "BH", alpha = alphaValue,contrast=c("condition","Primary Tumor","Solid Tissue Normal"))
#Nota
# si tengo mas grupos para comaprar res.2 =  y armo el otro contrast

#El resultado es una matriz con:
#1. promedio
#2. log2FoldChange  = log2(treated / untreated)
#3. pValue 
#4. adjusted pValue (es lo que usamos para rechazar H0)

res

#Nota sobre pValues = NA puede ser debido a 
#1. If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
#2. If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. 
#3. If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. 



#***************************************************************
#summary results
#***************************************************************
#describmos la informacion sobre los test que usamos
mcols(res)$description

summary(res)

#vemos cuantos dieron significativos
sum(res$padj < alphaValue, na.rm=TRUE)
#***************************************************************
#exploring results
#***************************************************************
#plotMA
#Muestra el Log2 fold change rspecto a al promedio de counts normalizado de todas las muestras
# En rojo si el adjusted pValue es <alpha
# Outliers se muestran como tiangulos ya sea que estan up o down
plotMA(res, ylim=c(-2,2))

#nos permite de manera interactiva al tocar los puntos en el plot que gen estamos tocando y guardarlos
#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

#Plot counts
#nos permite ver para una variable determinada los counts que hay en cada grupo
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#***************************************************************
#annotations
#***************************************************************
ens.id <- row.names(DF_RESULTS)
#remove ensemble verion
ens.id<- gsub('\\..+$', '', ens.id )
DF_FINAL_RESULT = cbind("ensembl_gene_id"= ens.id,as.data.frame(res))
COUNT_MATRIX = cbind("ensembl_gene_id"= ens.id,assay(dds))

library("biomaRt")
#http://www.ensembl.org/info/data/biomart/biomart_r_package.html
mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# Use BioMart to return gene names for a list of Ensembl IDs

head(listAttributes(ensembl))

gene.names <- getBM(filters= "ensembl_gene_id",
                    values = ens.id,
                    attributes= c("ensembl_gene_id","hgnc_symbol",'chromosome_name','start_position','end_position'),
                    mart= mart)



DF_RESULT <- merge(x = gene.names, 
                   y = DF_FINAL_RESULT, 
                   by.x="ensembl_gene_id",
                   by.y="ensembl_gene_id", sort = FALSE)


COUNT_MATRIX_RESULT <- merge(x = gene.names, 
                             y = COUNT_MATRIX, 
                             by.x="ensembl_gene_id",
                             by.y="ensembl_gene_id", sort = FALSE)

#***************************************************************
#Guardamos
#***************************************************************
#guardamos count matrix
write.csv(as.data.frame(COUNT_MATRIX_RESULT), file=paste0(HTCountsDetinyPath,"TCGA_COAD_CountMatrix.csv"))


#guardamos los resultados LFC original
write.csv(as.data.frame(DF_RESULT),file=paste0(HTCountsDetinyPath,"TCGA_COAD_LFC_results.csv"))


#Guardamos solo los que dieron significativo y ordenados por adjusted pValue
resSig <- subset(DF_RESULT, padj < alphaValue)
resOrdered <- resSig[order(resSig$padj),]
resOrdered
write.csv(resOrdered,file=paste0(HTCountsDetinyPath,"TCGA_COAD_LFC_SOLO_SIGNIFICATIVOS.csv"))

#***************************************************************
#Lista de enzimas
#***************************************************************
#filtramos la lista de enzimas y guardamos todos y por otro lado solo los significativos
#leemos la lista de enzimas
library(readr)
enzimas <- read_delim("../resources/enzimas.csv", 
                      ";", escape_double = FALSE, trim_ws = TRUE)

lista = cbind( Group = enzimas$Group,hgnc_symbol = toupper(enzimas$GeneSymbol))

filtradoEnzimas <- merge(x = lista  , 
                         y = DF_RESULT, 
                         by.x="hgnc_symbol",
                         by.y="hgnc_symbol",sort = FALSE )

write.csv(filtradoEnzimas,file=paste0(HTCountsDetinyPath,"TCGA_COAD_LFC_FILTRADO_LISTA_ENZIMAS.csv"))

#***************************************************************
#Reporting tools
#***************************************************************

des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',title = 'TCGA-COAD RNA-seq analysis of differential expression using DESeq2', reportDirectory = HTCountsDetinyPath)
publish(filtradoEnzimas,des2Report,  reportDir=HTCountsDetinyPath)
finish(des2Report)


library(GOstats)
goParams <- new("GOHyperGParams", geneIds = selectedIDs, universeGeneIds = universeIDs, annotation ="org.Mm.eg" , ontology = "MF",conditional = TRUE, testDirection = "over")
goResults <- hyperGTest(goParams)


goReport <- HTMLReport(shortName = 'TCGA-COAD GO ANALYSIS',title = "GO analysis of mockRnaSeqData",reportDirectory = HTCountsDetinyPath)
publish(goResults, goReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg.db")
finish(goReport)


library(Category)
params <- new("PFAMHyperGParams", geneIds= selectedIDs, universeGeneIds=universeIDs,annotation="org.Hs.eg.db",testDirection="over")
PFAMResults <- hyperGTest(params)

PFAMReport <- HTMLReport(shortName = 'TCGA-COAD PFAM ANALYSIS', title = "PFAM analysis of mockRnaSeqData", reportDirectory = HTCountsDetinyPath)
publish(PFAMResults, PFAMReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg.db",categorySize=5)
finish(PFAMReport)

indexPage <- HTMLReport(shortName = "indexRNASeq", title = "Analysis of TCGA-COAD",reportDirectory = HTCountsDetinyPath)
publish(Link(list(deReport,des2Report, goReport, PFAMReport), report = indexPage),
        indexPage)

