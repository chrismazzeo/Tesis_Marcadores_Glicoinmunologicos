#TODO Guardar modelo, error matrix y graficos
#TODO definir que quiero predecir tumor/ctrl  tumor stage y msi
#TODO ver si es con la vst y si es como hago con otros?
#TODO probarlo con el validation dataset
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


