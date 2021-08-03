list.of.packages <- c("readr","pheatmap","RColorBrewer","genefilter","gridExtra","mvShapiroTest")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#source("https://bioconductor.org/biocLite.R")

library(readr)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(gridExtra)
library(mvShapiroTest)
library(MVN)
#***************************************************************
#heatmaps
#***************************************************************
#***************************************************************
#heatmap 
#***************************************************************
#VST Gene list Expression accross samples
#vemos la varianza de cada gen para cada muestra - el promedio para ese gen
#si armamos clusters vemos  tb  coexpression de genes
#clustering de columnas muestra que samples se agrupan basados en la expresion d elos genes,uso como distancia "euclidean"
#clustering en las filas muestra coexpresion de genes, uso como  "correlation" o "euclidean"
#***************************************************************
samples <- read_delim(input.htSeqSampleFileName,"\t", escape_double = FALSE, trim_ws = TRUE)
vstMatrix <- read_csv(output.vstMatrix.FilteredByGeneLisAndSignificantFileName)


dim(vstMatrix)
vstMatrix[1:3,1:7]
vsdData = vstMatrix[,7:dim(vstMatrix)[2]]
vsdData = as.matrix(vsdData)
rownames(vsdData) = vstMatrix$EnsembleId 
colnames(vsdData) = samples$`File ID`
mat  <- vsdData - rowMeans(vsdData)


anno_col <- data.frame( Sample = samples$`Sample Type` )
rownames(anno_col) <-samples$`File ID`

anno_row <- data.frame( Group = vstMatrix$Group)
rownames(anno_row) <-vstMatrix$EnsembleId

group.color = list(Group = c(a=settings.color[1:length(levels(as.factor(vstMatrix$Group)))]), Sample = c("green", "pink"))
names(group.color$Group) = (levels(as.factor(vstMatrix$Group)))
names(group.color$Sample) = (levels(as.factor(samples$`Sample Type`)))

pheatmap(mat,
         main = paste0("VST-transformed values \n ",dim(vsdData)[1]," Genes Filtered by Glyco list & significant genes with adjusted pValue <",settings.de.alphaValue, "\n  rows = genes, clustering method = Correlation \n columns = Sample, clustering method = Euclidean\n "), cellwidth = settings.heatmap.cellWidth, cellheight = settings.heatmap.cellHeight *2,
         cluster_cols = TRUE, fontsize_col =  6, show_colnames = FALSE,
         cluster_rows = TRUE, fontsize_row = 6, 
         clustering_distance_columns = "euclidean",
         clustering_distance_rows = "correlation",
         annotation_col= anno_col,
         annotation_row = anno_row,
         annotation_colors = group.color,
         na_col= settings.heatmap.naColor,
         labels_row = vstMatrix$GeneSymbol,
         filename = output.heatmap.vst.GeneFilteredSifnificantExpressionFileName)
#***************************************************************
#correlation 
hacemos clustering metodo = correlacion de la count matrix y graficamos con heatmap, aca hay que hacer z score o log fold? ver como da
luego calculamos la correlacion sobre la count matrix, (no clusteamos sobre la correlacion porque no tiene sentido)aca podemos hacer heatmap o filtramos los que sean mayor <06 y pvalue significativo heatmap pero solo la diagonal inferior con algun package_version()
tb a modo de estar seguro podemos hacer lo msmo con la vst


#***************************************************************

a = read_csv(output.countMatrix.FileName)
a = as.data.frame(a)
dim(a)
a[1:3,1:4]
a = a[,-2:-3]
write_csv(a,output.infiltrationMatrixFileName)

countMatrix <- read_csv(output.countMatrix.FilteredByGeneListAnsSifnificantFileName)
dim(countMatrix)
countMatrix[1:3,1:7]
matrix = countMatrix[,7:dim(countMatrix)[2]]
matrix = as.matrix(matrix)
rownames(matrix) = countMatrix$EnsembleId
colnames(matrix) = samples$`File ID`
anno_col <- data.frame( Group = countMatrix$Group)
rownames(anno_col) <-countMatrix$EnsembleId


tmat = t(matrix)
shapiroResult = mvShapiro.Test(tmat)

if (shapiroResult$p.value>=0.05){
  correlationMethod ="Correlation method Pearson"
  correlation = cor(tmat)
  
}else{
  correlationMethod ="Correlation method Spearman"
  correlation = cor (tmat, method = "spearman")
}
sink(settings.correlation.normalityResultFileName)
print(shapiroResult$method)
print(shapiroResult$statistic)
print(shapiroResult$p.value)
print(correlationMethod)
sink(type = "message")
sink()

group.color = list(Group = c(a=settings.color[1:length(levels(as.factor(countMatrix$Group)))]))
names(group.color$Group) = (levels(as.factor(countMatrix$Group)))


pheatmap(correlation, 
         main = paste0("CountMatrix  Correlation  \n","(",dim(tmat)[2]," Genes Filtered by Glyco list & significant genes with adjusted pValue < ",settings.de.alphaValue,")"),
         cellwidth = settings.heatmap.cellWidth, cellheight = settings.heatmap.cellHeight *2,
         cluster_cols = FALSE, fontsize_col =  6, 
         cluster_rows = TRUE , fontsize_row = 6, 
         breaks = seq(from=-1, to=1, by=0.05),
         color = colorRampPalette(settings.heatmap.correlation)(40),
         clustering_distance_columns = "correlation",
         clustering_distance_rows = "correlation",
         annotation_row = anno_col,
         annotation_colors = group.color,
         labels_row = countMatrix$GeneSymbol,
         labels_col = countMatrix$GeneSymbol,
         filename = output.heatmap.countMatrix.CorrelationGeneFilteredCorrelationFileName)

correlation.cutoff= correlation
correlation.cutoff[abs(correlation.cutoff)<settings.correlation.cutoff] = 0
pheatmap(correlation.cutoff, 
         main = paste0("CountMatrix \n Correlation >",settings.correlation.cutoff,"\n (",dim(correlation.cutoff)[1]," Genes Filtered by Glyco list & significant genes with adjusted pValue < ",settings.de.alphaValue," )"),
         cellwidth = settings.heatmap.cellWidth, cellheight = settings.heatmap.cellHeight *2,
         cluster_cols = FALSE, fontsize_col =  6, 
         cluster_rows = TRUE , fontsize_row = 6, 
         breaks = seq(from=-1, to=1, by=0.05),
         color = colorRampPalette(settings.heatmap.correlation)(40),
         clustering_distance_columns = "correlation",
         clustering_distance_rows = "correlation",
         annotation_row = anno_col,
         annotation_colors = group.color,
         labels_row = countMatrix$GeneSymbol,
         labels_col = countMatrix$GeneSymbol,
         filename = output.heatmap.countMatrix.CorrelationGeneFilteredCorrelationHighCorrelationFileName)


 aca hay que filtrart la count matrix para los nombres de los genes porque esta mal
 guardar la matrix
 
correlation.high = correlation.cutoff[rowSums(correlation.cutoff) != 1,]
correlation.high = correlation.high[,colSums(correlation.high) != 0]


pheatmap(correlation.high, 
         main = paste0("Gene Filtered With High Correlation  \n ","(",dim(correlation.high)[1]," Gene Filtered by list & significant genes with adjusted pValue < ",settings.de.alphaValue," & high correlation > ",settings.correlation.cutoff,")"),
         cellwidth = settings.heatmap.cellWidth, cellheight = settings.heatmap.cellHeight *2,
         cluster_cols = TRUE, fontsize_col =  6, show_colnames = FALSE,
         cluster_rows = TRUE, fontsize_row = 6, 
         breaks = seq(from=-1, to=1, by=0.05),
         color = colorRampPalette(settings.heatmap.correlation)(40),
         filename = output.heatmap.countMatrix.CorrelationGeneFilteredCorrelationHighCorrelationSummaryFileName)

uniPlot(correlation.high, type="histogram")
#******







































#1 heatmap DE result todos,  los padj> o na(padjs)  van como lfc = 0
#***************************************************************
DE_ResultMatrix.heatmap <- read_csv(output.DE.FileName)
DE_ResultMatrix.heatmap$log2FoldChange[DE_ResultMatrix.heatmap$padj>=settings.de.alphaValue | is.na(DE_ResultMatrix.heatmap$padj)  ] = 0
DE_ResultMatrix.heatmap$log2FoldChange[abs(DE_ResultMatrix.heatmap$log2FoldChange)<1] = 0

#convertimos una columna en una matrix
lfc = DE_ResultMatrix.heatmap$log2FoldChange
height = ceiling(length(lfc)/settings.heatmap.maxColumn)
matrix = matrix(NA, nrow = height, ncol=settings.heatmap.maxColumn)
dim(matrix)
for (i in 1:height){
  for (j in 1:settings.heatmap.maxColumn){
    matrix[i,j] =  lfc[(i-1)*settings.heatmap.maxColumn + j]
    
  }  
}

pheatmap(matrix,  
         main = "Differential expression results",
         cluster_rows = FALSE, cluster_cols = FALSE,cellwidth = 7,cellheight = 2,
         na_col= "black",color = c("red",),
         filename = output.heatmap.deResultsFileName)

hay que poner las legenda cuantos up y cuantos down
Separo por tipo de gen?????

# 2 opciones https://support.bioconductor.org/p/64146/
# si hago heatmap de shrunken LFC no hay problema con los genes que no dieron significativo en adjPvalue 
#si hago directo el heatmp del resultado de deseq2, los que no dieron significados tengo que asiganarle NA, o filtrar y quedarme con los que dieron significativos
DE_ResultMatrix.heatmap <- read_csv(output.DE.FileName)
dim(DE_ResultMatrix.heatmap)
DE_ResultMatrix.heatmap = DE_ResultMatrix.heatmap[!is.na(DE_ResultMatrix.heatmap$padj),]
dim(DE_ResultMatrix.heatmap)
DE_ResultMatrix.heatmap$log2FoldChange[DE_ResultMatrix.heatmap$padj>=alphaValue] = NA
#DE_ResultMatrix.heatmap = DE_ResultMatrix.heatmap[order(DE_ResultMatrix.heatmap$log2FoldChange, decreasing = TRUE),]


lfc = DE_ResultMatrix.heatmap$log2FoldChange
lfc = na.omit(lfc)

height = ceiling(length(lfc)/heatmap.width)

DE_Matrix = matrix(NA, nrow = height, ncol=heatmap.width)
dim(DE_Matrix)
for (i in 1:height){
  for (j in 1:heatmap.width){
    DE_Matrix[i,j] =  lfc[(i-1)*heatmap.width + j]
 
  }  
}
sum(is.na(DE_Matrix))

pheatmap(DE_Matrix,  
         main = paste0("\nDifferential expression \n(",length(lfc)," significant genes with adjusted pValue < ",alphaValue,")\n"),
         cluster_rows = FALSE, cluster_cols = FALSE,cellwidth = heatmap.cellWidth,cellheight = heatmap.cellHeight,
         color = colorRampPalette(heatmap.colorRamp)(max(lfc)-min(lfc)),na_col= heatmap.naColor,
         filename = paste0(heatmap.path,heatmap.DE_RESULTS))
#***************************************************************
#heatmap DE result significant lfc> lfcThreshold
#***************************************************************
DE_ResultMatrix.heatmap <- read_csv(paste0(folder.results.alpha.DE,output.deLFCCutOffFileName))
dim(DE_ResultMatrix.heatmap)
DE_ResultMatrix.heatmap = DE_ResultMatrix.heatmap[!is.na(DE_ResultMatrix.heatmap$padj),]
dim(DE_ResultMatrix.heatmap)
DE_ResultMatrix.heatmap$log2FoldChange[DE_ResultMatrix.heatmap$padj>=alphaValue] = NA
#DE_ResultMatrix.heatmap = DE_ResultMatrix.heatmap[order(DE_ResultMatrix.heatmap$log2FoldChange, decreasing = TRUE),]


lfc = DE_ResultMatrix.heatmap$log2FoldChange
lfc = na.omit(lfc)
lfc = lfc[abs(lfc)>lfcThreshold]

height = ceiling(length(lfc)/heatmap.width)

DE_Matrix = matrix(NA, nrow = height, ncol=heatmap.width)
dim(DE_Matrix)
for (i in 1:height){
  for (j in 1:heatmap.width){
    DE_Matrix[i,j] =  lfc[(i-1)*heatmap.width + j]
    
  }  
}
sum(is.na(DE_Matrix))

pheatmap(DE_Matrix,  
         main = paste0("\nDifferential expression \n(",length(lfc)," significant genes with adjusted pValue < ",alphaValue," & |LFC| > ", lfcThreshold,")\n"),
         cluster_rows = FALSE, cluster_cols = FALSE,cellwidth = heatmap.cellWidth,cellheight = heatmap.cellHeight,
         color = colorRampPalette(heatmap.colorRamp)(max(lfc)-min(lfc)),na_col= heatmap.naColor,
         filename = paste0(heatmap.path,heatmap.DE_RESULTS_LFC))
#***************************************************************
#heatmap DE result filtered by gene list 
#***************************************************************
DE_ResultMatrix.heatmap <- read_csv(paste0(folder.results.alpha.DE,output.deGeneFilteredFileName))
dim(DE_ResultMatrix.heatmap)
#DE_ResultMatrix.heatmap = DE_ResultMatrix.heatmap[!is.na(DE_ResultMatrix.heatmap$padj),]
#dim(DE_ResultMatrix.heatmap)
#DE_ResultMatrix.heatmap$log2FoldChange[DE_ResultMatrix.heatmap$padj>=alphaValue] = NA
#DE_ResultMatrix.heatmap = DE_ResultMatrix.heatmap[order(DE_ResultMatrix.heatmap$log2FoldChange, decreasing = TRUE),]


lfc = DE_ResultMatrix.heatmap$log2FoldChange
heatmap.width2 = 20
height = ceiling(length(lfc)/heatmap.width2)

DE_Matrix = matrix(NA, nrow = height, ncol=heatmap.width2)
dim(DE_Matrix)
for (i in 1:height){
  for (j in 1:heatmap.width2){
    DE_Matrix[i,j] =  lfc[(i-1)*heatmap.width2 + j]
    
  }  
}
sum(is.na(DE_Matrix))
pheatmap(DE_Matrix,  
         main = paste0("\nDifferential expression \n(",length(lfc)," significant genes with adjusted pValue < ",alphaValue,", Filtered By gene list)\n"),fontsize = 4,
         cluster_rows = FALSE, cluster_cols = FALSE,cellwidth = 20,cellheight = 7,
         color = colorRampPalette(heatmap.colorRamp)(max(lfc)-min(lfc)),na_col= heatmap.naColor,
         filename = paste0(heatmap.path,heatmap.DE_RESULTS_Filtered))

#nx1
pheatmap(DE_ResultMatrix.heatmap$log2FoldChange,  
         fontsize = 4,
         cluster_rows = FALSE, cluster_cols = FALSE,cellwidth =heatmap.width2,cellheight = 7,
         show_colnames = TRUE, show_rownames = TRUE,
         labels_row =  DE_ResultMatrix.heatmap$GeneSymbol,
         color = colorRampPalette(heatmap.colorRamp)(max(lfc)-min(lfc)),na_col= heatmap.naColor,
         filename =paste0(heatmap.DE_RESULTS_FilteredListxGroup,i)) 
#***************************************************************
#heatmap DE result filtered by gene list x family group
#***************************************************************

familys = levels(as.factor(DE_ResultMatrix.heatmap$Group))
plot_list=list()
for (i in familys){
   group = subset(DE_ResultMatrix.heatmap,DE_ResultMatrix.heatmap$Group == i)
   groupMatrix = matrix(group$log2FoldChange)
  if (dim(groupMatrix)[1]>=2){
    print (i)
    pheatmap(groupMatrix,  
             fontsize = 8,fontsize_row = 6,
             main = i,silent = FALSE,
             cluster_rows = FALSE, cluster_cols = FALSE,cellwidth = 20,cellheight = 7,
             show_colnames = FALSE, show_rownames = TRUE,
             labels_col = i,labels_row =  group$GeneSymbol,
             color = colorRampPalette(heatmap.colorRamp)(max(DE_ResultMatrix.heatmap$log2FoldChange)-min(DE_ResultMatrix.heatmap$log2FoldChange)),na_col= heatmap.naColor,
             filename = paste0(heatmap.DE_RESULTS_FilteredListxGroup,gsub("/", " ", i),".pdf"))
     
  }
}  



