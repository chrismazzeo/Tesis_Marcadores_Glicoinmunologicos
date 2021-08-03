source("./src/shared/packages.r")
#*********************************************************************************************************
project.resultsPath = "./results/ConsolidadoAll/"
geneListFilterPathHS = "./resources/GlycoGeneList/Glyco_HS_GeneList.csv"
geneListFilterPathMM = "./resources/GlycoGeneList/Glyco_MM_GeneList.csv"
geneListHS =   HELPER_LoadCSV(geneListFilterPathHS)
geneListMM =   HELPER_LoadCSV(geneListFilterPathMM)
tcgaResultPath = "./Results/human/gdc/TCGA/"
mouseResultPath = "./Results/mouse/"
dir.create(project.resultsPath, recursive = TRUE )
#*********************************************************************************************************
#heatmap settings

heatmap.clearEmptyRows = TRUE
heatmap.breakUp = 10
heatmap.breakDown = 10
heatmap.rowFontSize = 6
heatmap.colFontSize = 6
heatmap.celWidth = 20
heatmap.cellHeight = 10
heatmap.colorUp = colorRampPalette((brewer.pal(n = 8, name = "Reds")))(heatmap.breakUp)
heatmap.colorDown = colorRampPalette(rev(brewer.pal(n = 8, name = "Blues")))(heatmap.breakDown)
heatmap.colorZero = "white"
heatmap.colorNa = "white"
heatmap.colors=c( heatmap.colorDown,heatmap.colorZero,heatmap.colorUp)

#*********************************************************************************************************
#Hay que poder pasarle para que lea el que es todo significativos y corregir donde dice ensemblID
#nHay que agregar el genename para colaparslo

tcgaRes = readAndMergeFiles (resultPath = tcgaResultPath,
                             geneList = geneListHS ,
                             outputDir = project.resultsPath,
                             filePrefix = "HS")


mouseRes = readAndMergeFiles (resultPath = mouseResultPath,
                              geneList = geneListMM,
                              outputDir = project.resultsPath,
                              filePrefix = "MM",
                              mergeByGeneName = TRUE)



#collapse into genes
tcgaRes.logfc0 = collapseGenesDataFrame(tcgaRes$H0log2FC0)
tcgaRes.logfc1 = collapseGenesDataFrame(tcgaRes$H0log2FC1)

mouseRes$H0log2FC0  = mouseRes$H0log2FC0[,-3]
mouseRes$H0log2FC1  = mouseRes$H0log2FC1[,-3]

mouseRes.logfc0 = collapseGenesDataFrame(mouseRes$H0log2FC0 )
mouseRes.logfc1 = collapseGenesDataFrame(mouseRes$H0log2FC1)
mouseRes.logfc0$geneName = toupper(mouseRes.logfc0$geneName )
mouseRes.logfc1$geneName = toupper(mouseRes.logfc1$geneName )


#Log2FC 0 

finalLogFc0 = merge(mouseRes.logfc0, tcgaRes.logfc0, by.x = "geneName", by.y = "geneName", all = TRUE )
rownames(finalLogFc0) = finalLogFc0$geneName
extraData = data.frame(Genes = finalLogFc0$geneName, Family = finalLogFc0$family.x)

finalLogFc0 = finalLogFc0[,-c(1,22, 81)]
for (i in 1:dim(finalLogFc0)[2]){
  finalLogFc0[,i] = as.numeric(finalLogFc0[,i] ) 
}
finalLogFc0 = rapply( finalLogFc0, f=function(x) ifelse(is.na(x),0,x), how="replace" )

annotationCol = data.frame(Type = rep("Mouse",dim(mouseRes.logfc0)[2]-2))
annotationCol= rbind(annotationCol, data.frame(Type = rep("TCGA",dim(tcgaRes.logfc0)[2]-2)))
rownames(annotationCol) = colnames(finalLogFc0)
annotationRow =data.frame(Family = extraData$Family)
rownames(annotationRow) = extraData$Genes

#TODO cambiar nombres a las comapraciones y orden
#TODO filtrar columnas en 0 y genes y familias?
finalMerge = cbind(extraData,finalLogFc0)
HELPER_SAVE_DATA_FRAME(finalMerge, paste0(project.resultsPath,"FinalMerge.csv"))


heatmap.breaks = unique(c(
  seq(min(finalLogFc0, na.rm = TRUE),-1, length=heatmap.breakDown),
  0,
  rev(seq(max(finalLogFc0, na.rm = TRUE),1, length=heatmap.breakUp))
))


pheatmap(finalLogFc0, 
         cellwidth = heatmap.celWidth, cellheight =heatmap.cellHeight,fontsize_number = heatmap.rowFontSize,fontsize_row = heatmap.colFontSize,
         na_col =heatmap.colorNa, color=heatmap.colors,breaks = heatmap.breaks, 
         annotation_row = annotationRow, annotation_col = annotationCol,
         
         display_numbers = FALSE, 
         cluster_rows = FALSE, cluster_cols = FALSE,  show_rownames = TRUE, show_colnames = TRUE,
         filename = paste0(project.resultsPath,"finalLog2FC0.pdf"))

cor = HELPER_RowNamesAsFirstColumn(correlationResult$r, "Comparison")
padj = HELPER_RowNamesAsFirstColumn(correlationResult$p, "Comparison")
HELPER_SAVE_DATA_FRAME(cor,paste0(project.resultsPath,"correlation_R.csv"))
HELPER_SAVE_DATA_FRAME(padj,paste0(project.resultsPath,"correlation_pValue.csv"))

#correlation
correlationResult = corr.test(x = finalLogFc0,use = "pairwise", method = "spearman", adjust = "BH", alpha=0.05, ci = FALSE)
pdf(paste0(project.resultsPath,"correlation.pdf"), width = 30, height= 30)
corrplot(corr = correlationResult$r,
         col = heatmap.colors,
         title ="Correlation", 
         diag = FALSE,
         type = "lower",
         tl.col="black", #Text label color and rotation
         number.cex =0.4,
         cl.cex = 1,
         p.mat = correlationResult$p, sig.level = 0.05,insig = "blank",
         na.label = "  ",
         mar=c(0,0,1,0)
) 

dev.off()


#verificar los valores para ver que no meti la pata con tantas operaciones

#remove empty rows and cols
# if (heatmap.clearEmptyRows){
#   dim(collapsedMatrix.matrix)
#   removeCols = colSums(collapsedMatrix.matrix, na.rm =TRUE) == 0
#   removeRows = rowSums(collapsedMatrix.matrix, na.rm =TRUE) == 0
#   
#   collapsedMatrix.matrix = collapsedMatrix.matrix[!removeRows,!removeCols]
#   
#   annotation = annotation[!removeRows,]
#   annotation = droplevels(annotation)
#   
#   geneName = geneName[!removeRows]
#   geneName = droplevels(geneName)
#   dim(collapsedMatrix.matrix)
# }

#*******************************************************************************************************
#feature selection y predictions
#*******************************************************************************************************

#https://www.biostars.org/p/86981/
source("./src/shared/packages.r")
source("./src/shared/settings.r")
library(mltest)
library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)
require(reshape2)
library(pROC)
#*******************************************************************************************************
useVSD = TRUE
filterGeneList = TRUE
genesTop = 5
outputDir= "./results/consolidado/predicciones/"
dir.create(outputDir)
dir.create(file.path( outputDir), recursive = TRUE, showWarnings = FALSE)

project.resourcePath = "./resources/human/GDC/TCGA/COAD/" 
vsdPath = "./Results/Human/GDC/TCGA/COAD_ALL_VARS/HTSEQ/CountMatrix/vsd.csv"
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList_filtered.csv"
geneListFilterPath = "./resources/GlycoGeneList/Glyco_HS_GeneList.csv"

samplesList = HELPER_LoadCSV(samplesListPath)
geneList =   HELPER_LoadCSV(geneListFilterPath)
vsd = HELPER_LoadTSV(vsdPath)
vsd = vsd[,-1]

varName = "GROUP"
#group = data.frame(Labels = paste(samplesList$tumor_stage ,samplesList$msi_status, sep ="_"), stringsAsFactors = FALSE)
group = paste(samplesList[,varName], sep ="_")
groups = gsub(" ","_",group)
groups = gsub("-","_",groups)
groups[samplesList$GROUP == "Solid Tissue Normal"] = "Control"
groups = as.factor(groups)
groups

#*******************************************************************************************************

#samplesList = samplesList[samplesList$GROUP == "Primary Tumor",]

htseq.rawCountMatrix = DESeq2_Merge_TCGA_ZippedExpressionFiles(sourceFolder = paste0(project.resourcePath,"rawData/htseq"),
                                                               sampleFolder = samplesList$File_ID,
                                                               samplesFileName =samplesList$File_Name,
                                                               sampleID = samplesList$Sample_ID,
                                                               outputPath = NULL)
rownames(vsd) = rownames(htseq.rawCountMatrix)
colnames(vsd) = colnames(htseq.rawCountMatrix)
dim(htseq.rawCountMatrix)
vsd = t(vsd)
dim(vsd)

htseq.rawCountMatrix = t(htseq.rawCountMatrix)
dim(htseq.rawCountMatrix)


#*******************************************************************************************************
#filtramos por la lista de glyco y elgirmos vsd o htseq
#*******************************************************************************************************
if (filterGeneList){
  htseq.rawCountMatrix.filtered = htseq.rawCountMatrix[,colnames(htseq.rawCountMatrix) %in% geneList$ensembl_gene_id]
  vsd.filtered = vsd[,colnames(vsd) %in% geneList$ensembl_gene_id]
}else{
  htseq.rawCountMatrix.filtered = htseq.rawCountMatrix
  vsd.filtered = vsd
}

if (useVSD){
  matrix = vsd.filtered 
}else{
  matrix = htseq.rawCountMatrix.filtered
}

#*******************************************************************************************************
#Filtramos las variables
#Solos las que tengan <.7 coeficiente de variacion para reducir la cantidad de variables
#*******************************************************************************************************
# dim(matrix)
# table(colcvs(matrix)>0.2)
# matrix = matrix[,colcvs(matrix)>0.2]
# table(colcvs(matrix))
#*******************************************************************************************************
# Random Forest
#*******************************************************************************************************
matrix = matrix[groups != "NA",]
groups = groups[groups != "NA"]
groups = droplevels(groups)

rf_output=randomForest(x=matrix, y=groups, importance = TRUE, ntree = 10001, proximity=TRUE,  na.action = na.omit)


save(rf_output, file= paste0(outputDir,"RF_model_",varName))


#*******************************************************************************************************
#peso de cada variable
#ploteamos a ver como sale los primeros 5 genes
#*******************************************************************************************************
rf_importances=importance(rf_output, scale=FALSE)
rf_importances = as.data.frame(rf_importances)
rf_importances = rf_importances[order(-rf_importances$MeanDecreaseGini ),] 
#View(rf_importances)

selectedGenesID = rownames(rf_importances)[1:genesTop]
selectedGenesExpression = matrix[,colnames(matrix) %in% selectedGenesID]
selectedGenesExpression = data.frame(Group = groups, selectedGenesExpression)


selectedGenesExpression.melt = melt(selectedGenesExpression ,id.vars = "Group")

#*******************************************************************************************************
#Error matrix
#*******************************************************************************************************

confusion=rf_output$confusion
confusion
confusion = HELPER_RowNamesAsFirstColumn(confusion,"Group")
HELPER_SAVE_DATA_FRAME(confusion, paste0(outputDir, "ErrorMatrix_",varName,".csv"))

errorMatrix = ml_test(rf_output$predicted, rf_output$y, output.as.table = TRUE)
errorMatrix = HELPER_RowNamesAsFirstColumn(errorMatrix,"Group")
HELPER_SAVE_DATA_FRAME(errorMatrix, paste0(outputDir, "ErrorMatrixStats_",varName,".csv"))
#*******************************************************************************************************
# AUC ROC
#*******************************************************************************************************
predictions = as.numeric(predict(rf_output))
roc.multi = multiclass.roc(groups, predictions, percent=TRUE)

rs <- roc.multi[['rocs']]

pdf(paste0(outputDir, varName,".pdf"))
ggplot(data = selectedGenesExpression.melt, aes(x=variable, y=value)) + geom_boxplot(aes(fill=Group))+
  settings.graphics.theme +  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  scale_fill_manual(values=c(settings.graphics.colors[1:length(levels(groups))]))

plot(selectedGenesExpression, col = selectedGenesExpression$Group, main = varName)
varImpPlot(rf_output, type=2,sort = TRUE ,n.var=genesTop, scale=FALSE, main=paste0("Variable Importance (Gini) for top ",genesTop," predictors"))
MDSplot(rf_output, groups, xlab="", ylab="", palette= settings.graphics.colors[1:length(levels(groups))], main="MDS plot")

plot.roc(rs[[1]])
sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))
comp =sapply(1:length(rs),function(i) paste0(rs[[i]]$levels[1],"/",rs[[i]]$levels[2]))
legend("bottomright", legend=comp,lwd=2, col =seq(1:length(a)))
dev.off()

# out <- histbackback(split(rf_output$votes, groups), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for patients classified by RF', axes=TRUE, ylab="Fraction votes")
# barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE)
# barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE)

#load("RF_model_Type")
#TODO predecir el resto del dataframe






