list.of.packages <- c("readr","dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#source("https://bioconductor.org/biocLite.R")

library(readr)
library(dplyr)



#***************************************************************
#Load 
#***************************************************************

countMatrix <- read_csv(paste0(folder.results, input.countMatrixFileName))
countMatrix = countMatrix[,-1]
countMatrix[1:3,1:3]

vstMatrix =read_csv(paste0(folder.results, output.vstMatrixSortedBySampleType))
vstMatrix = countMatrix[,-1]
vstMatrix[1:3,1:3]

sampleTypes = read_delim(paste0(folder.resources, input.htSeqSampleFileName),"\t", escape_double = FALSE, trim_ws = TRUE)
sampleTypes[1:3,1:3]

DEMatrix <- read_csv(paste0(folder.results, input.deFileName))
DEMatrix = DEMatrix[,-1]
DEMatrix[1:3,1:3]

enzimas <- read_delim(input.geneFilteringFileName, 
                      ";", escape_double = FALSE, trim_ws = TRUE)
lista = cbind( Group = enzimas$Group,hgnc_symbol = toupper(enzimas$GeneSymbol))


#***************************************************************
#sort samples
#***************************************************************
sampleTypes.sorted =sampleTypes[order(sampleTypes$`Sample Type`,  decreasing = TRUE), ]
fileNames = strsplit(sampleTypes.sorted$`File Name`, ".gz")
sampleTypes.sorted$`File Name` = unlist(fileNames)
write.csv(sampleTypes.sorted,file=paste0(folder.results, output.htSeqSampleSortedFileName))

#***************************************************************
# Differential expression 
#***************************************************************
#***************************************************************
# Solo SIGNIFICATIVOS 
#***************************************************************
dim(DEMatrix)
DESignificativos <- subset(DEMatrix, padj < alphaValue)
dim(DESignificativos)
DESignificativos[1:3,1:6]
write.csv(DESignificativos,file=paste0(folder.results, output.deSignificantFileName))
#***************************************************************
#log2 CUTOFF de significativos 50% up o down
#***************************************************************
DESignificativosCutOff <- subset(DESignificativos, DESignificativos$log2FoldChange<=-logFoldChangeCutOff | DESignificativos$log2FoldChange>=logFoldChangeCutOff)
dim(DESignificativosCutOff)
write.csv(DESignificativosCutOff,file=paste0(folder.results, output.deLFCCutOffFileName))
#***************************************************************
#Filtrado x Lista de enzimas
#***************************************************************
filtradoEnzimas <- merge(x = lista  , 
                         y = DEMatrix, 
                         by.x="hgnc_symbol",
                         by.y="hgnc_symbol",sort = FALSE )

write.csv(filtradoEnzimas,file=paste0(folder.results,output.deGeneFilteredFileName))

#***************************************************************
#VST Matrix
#***************************************************************
#***************************************************************
#VST todas las muestras Sorted by sample type
#***************************************************************
vstMatrix.sorted = as.data.frame(vstMatrix)
vstMatrix.sorted = vstMatrix[,-1:-4]
vstMatrix.sorted = vstMatrix.sorted[,order(match(names(vstMatrix.sorted),paste0(sampleTypes.sorted$`Case ID`,".",sampleTypes.sorted$`File Name`)), decreasing = FALSE)]
vstMatrix.sorted[1:3,1:10]
dim(vstMatrix.sorted)
write.csv(vstMatrix.sorted,file=paste0(folder.results,output.vstMatrixSortedBySampleType))
#***************************************************************
# todas las muestras Filtrado x Lista de enzimas 
#***************************************************************
vstMatrix.sorted = cbind(hgnc_symbol = vstMatrix$hgnc_symbol,vstMatrix.sorted)
filtradoEnzimas <- merge(x = lista  , 
                         y = vstMatrix.sorted, 
                         by.x="hgnc_symbol",
                         by.y="hgnc_symbol",sort = FALSE )

write.csv(filtradoEnzimas,file=paste0(folder.results,output.vstMatrixGeneFiltered))


#***************************************************************
#Count Matrix
#***************************************************************
#***************************************************************
#Count MAtrix todas las muestras Sorted by sample type
#***************************************************************
countMatrix.sorted = as.data.frame(countMatrix)
countMatrix.sorted = countMatrix[,-1:-5]
countMatrix.sorted = countMatrix.sorted[,order(match(names(countMatrix.sorted),paste0(sampleTypes.sorted$`Case ID`,".",sampleTypes.sorted$`File Name`)), decreasing = FALSE)]
countMatrix.sorted[1:3,1:10]
dim(countMatrix.sorted)
write.csv(countMatrix.sorted,file=paste0(folder.results,output.countMatrixSortedBySampleType))
#***************************************************************
# todas las muestras Filtrado x Lista de enzimas 
#***************************************************************
countMatrix.sorted = cbind(hgnc_symbol = countMatrix$hgnc_symbol,countMatrix.sorted)
filtradoEnzimas <- merge(x = lista  , 
                         y = countMatrix.sorted, 
                         by.x="hgnc_symbol",
                         by.y="hgnc_symbol",sort = FALSE )

write.csv(filtradoEnzimas,file=paste0(folder.results,output.countMatrixGeneFiltered))






