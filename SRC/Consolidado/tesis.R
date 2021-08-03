####### WARNING ###################################
#boxplot with pvalues, hay que editar a mano y reemplazar por los padjusted
#################################################
library(ggsci)
library(RColorBrewer)
library(readr)
library(cowplot)
library(ggpubr)
library(pheatmap)
library(writexl)
library(psych)
library(corrplot)
library(stringr)
source("./src/shared/functions/HELPER_Function.r")
# install.packages("devtools")
# devtools::install_github("rlbarter/superheat")
# library(superheat)

inflamacion = HELPER_LoadCSV("./results/tesis/ratón/RatonFinalMergeTesis.csv")
colnames(inflamacion) = c("Genes",	
                     "Family",
                     "GSE43338: AOM/DSS VS Control",
                     "GSE43338: Apc min/+ (spCRC) VS Control",
                     "GSE57533: AOM/DSS VS Control",
                     "GSE57533: DSS VS Control",
                     "GSE64658: Late Distal AOM/DSS  VS Control"
                     )

inflamacion = inflamacion[,-4]
inflamacion[,3:6][abs(inflamacion[,3:6])<1.5] = NA
inflamacion = na.omit(inflamacion)
HELPER_SAVE_DATA_FRAME(inflamacion,"./results/tesis/ratón/inflamacion.csv")
rownames(inflamacion) = inflamacion$Genes
inflamacion = inflamacion[order(inflamacion[,3],decreasing = TRUE),]
pheatmap(inflamacion[,3:6], main = "Inflammation", 
         cellheight = 10,cellwidth = 20,
         cluster_cols = FALSE, cluster_rows = FALSE,
         treeheight_row = 0, treeheight_col = 0,
         fontsize_row = 7,
         filename = "./results/tesis/ratón/inflamacion.pdf")
# pdf("caca.pdf")
# 
# superheat(inflamacion[,3:6],
#           heat.pal = c( "blue","white","red"),heat.pal.values = c(0, 0.4,1),
#           bottom.label.text.angle = 90)
# dev.off()

cacrc = HELPER_LoadCSV("./results/tesis/ratón/RatonFinalMergeTesis.csv")
colnames(cacrc) = c("Genes",	
                          "Family",
                          "GSE43338: AOM/DSS VS Control",
                          "GSE43338: Apc min/+ (spCRC) VS Control",
                          "GSE57533: AOM/DSS VS Control",
                          "GSE57533: DSS VS Control",
                          "GSE64658: Late Distal AOM/DSS  VS Control"
)

cacrc = cacrc[,-c(3,5)]
cacrc[,3:5][abs(cacrc[,3:5])<1.5] = NA
cacrc = na.omit(cacrc)
HELPER_SAVE_DATA_FRAME(cacrc,"./results/tesis/ratón/cacrc.csv")
rownames(cacrc) = cacrc$Genes
cacrc = cacrc[order(cacrc[,3],decreasing = TRUE),]
pheatmap(cacrc[,3:5], main = "CACRC",
         cellheight = 10,cellwidth = 20,
         treeheight_row = 0, treeheight_col = 0, 
         cluster_cols = FALSE,cluster_rows = FALSE,
         fontsize_row = 7,
         filename = "./results/tesis/ratón/CACRC.pdf")


crc = HELPER_LoadCSV("./results/tesis/ratón/RatonFinalMergeTesis.csv")
colnames(crc) = c("Genes",	
                    "Family",
                    "GSE43338: AOM/DSS VS Control",
                    "GSE43338: Apc min/+ (spCRC) VS Control",
                    "GSE57533: AOM/DSS VS Control",
                    "GSE57533: DSS VS Control",
                    "GSE64658: Late Distal AOM/DSS  VS Control"
)
crc = crc[,-6]
crc[,3:6][abs(crc[,3:6])<1.5] = NA
crc = na.omit(crc)
HELPER_SAVE_DATA_FRAME(crc,"./results/tesis/ratón/crc.csv")
rownames(crc) = crc$Genes
crc = crc[order(crc[,3],decreasing = TRUE),]
pheatmap(crc[,3:6], main = "CRC", 
         cellheight = 10,cellwidth = 20,
         cluster_cols = FALSE,cluster_rows = FALSE,
         treeheight_row = 0, treeheight_col = 0, 
         fontsize_row = 7,
          filename = "./results/tesis/ratón/CRC.pdf")

#*************************************************************************************************************************************
#***OK REVISADO **********************************************************************************************************************
#*************************************************************************************************************************************
#immuCC GSE64658   HAy otra version 
#*************************************************************************************************************************************


project.immucc.svrPath =  "./resources/mouse/GSE64658/immucc/GSE64658.SVR.csv"
samplesListPath = "./resources/mouse/GSE64658/GSE64658_SampleList.csv"
samplesList = HELPER_LoadCSV(samplesListPath)
GSE64658.immucc =HELPER_LoadCSV(project.immucc.svrPath)
rownames(GSE64658.immucc)  = GSE64658.immucc$X
GSE64658.immucc = GSE64658.immucc[,-1]  
project.groupsToCompare = list( c("Distal_Colon_AOM_DSS_Early","Distal_Colon_Control"),
                             c("Distal_Colon_AOM_DSS_Late","Distal_Colon_Control"),
                             c("Proximal_Colon_AOM_DSS_Early","Proximal_Colon_Control"),
                             c("Proximal_Colon_AOM_DSS_Late","Proximal_Colon_Control"),
                             c("Distal_Colon_AOM_DSS_Late","Distal_Colon_AOM_DSS_Early"),
                             c("Proximal_Colon_AOM_DSS_Late","Proximal_Colon_AOM_DSS_Early"))



  
immuneMatrix =data.frame(GSE64658.immucc, Factor = samplesList$GROUP)
vars = colnames(immuneMatrix)[1:(dim(immuneMatrix)[2]-1)]
a = list()
testResults = data.frame()
for  (i in 1:length(vars)){

  data.formula = formula(paste(vars[i] ,"~ Factor"))
  testResults.temp =  compare_means(data.formula,comparisons = project.groupsToCompare, p.adjust.method = "fdr", method='wilcox.test',data = immuneMatrix)
  testResults = rbind(testResults,testResults.temp)
  
 #  pairwise.wilcox.test(x = immuneMatrix[,i], g =  immuneMatrix$Factor, p.adjust.method = "fdr")

  #esto muestra pvalue y no padjusted es una cagada
  p2 = ggboxplot(immuneMatrix,  x ="Factor",y = vars[i],
                 title = gsub("\\."," ",vars[i]),
                 subtitle = " ",
                 fill = "Factor", add = c("jitter"),width=0.5,
                 order = c( "Proximal_Colon_Control","Proximal_Colon_AOM_DSS_Early", "Proximal_Colon_AOM_DSS_Late",
                           "Distal_Colon_Control", "Distal_Colon_AOM_DSS_Early", "Distal_Colon_AOM_DSS_Late"),
                 xlab = "",ylab ="TIL (%)",
                 repel = TRUE,combine = FALSE)+
    

  stat_compare_means(comparisons = project.groupsToCompare,
                     method = "wilcox.test",
                     p.adjust.methods = "fdr"
                    ) 
  # +  stat_compare_means(label.y.npc = 1, label.x.npc = 0.3)  Kruskall

  #para agregar a mano padjusted, el tema es que te agrga todas las comapraciones y no las que queres  
  # stat_pvalue_manual(testResults, label = "p.adj",
  #                    y.position = seq(max(immuneMatrix[,i])+1, 10, length.out = dim(testResults)[1]))

  
  
   
 p2 = ggpar(p2, 
         font.title = c(14,"bold", "black"),palette = "npg"
        )
    p3 = p2
   p2 =  p2+rremove("x.text") +rremove("legend") +rremove("x.ticks")
   
   a[[i]] = p2
  #print(p2)  
 
  
}

HELPER_SAVE_DATA_FRAME(as.data.frame(testResults),"./results/tesis/ratón/immuCC/GSE64658/testResults.csv" )


p3  = p3 +rremove("x.title") +rremove("x.axis")+rremove("y.text") +rremove("y.ticks")+rremove("y.axis")+rremove("y.title") 
p3$labels$title = ""

p3 = ggpar(p3,
           legend = "right", legend.title = "Group", font.legend =c(12,"bold", "black"))
  a[[i+1]] = p3
  
b =   plot_grid(plotlist =a,
          ncol = 5, nrow = 6)

ggsave(b,file=paste0("./results/tesis/ratón/immuCC/GSE64658/SVR_BoxplotWithPValues.pdf"), width = 15, height =25, limitsize = FALSE)


#*************************************************************************************************************************************
#***OK REVISADO **********************************************************************************************************************
#*************************************************************************************************************************************
#immuCC GSE57533
#*************************************************************************************************************************************
project.immucc.svrPath =  "./resources/mouse/GSE57533/immucc/GSE57533.SVR.csv"
samplesListPath = "./resources/mouse/GSE57533/GSE57533_sampleList.csv"
samplesList = HELPER_LoadCSV(samplesListPath)
GSE57533.immucc =HELPER_LoadCSV(project.immucc.svrPath)
rownames(GSE57533.immucc)  = GSE57533.immucc$X
GSE57533.immucc = GSE57533.immucc[,-1]  
project.groupsToCompare.inmune = project.groupsToCompare  #levantrlo de GSE57533_setup.r
project.groupsToCompare = list( c("DSS","CONTROL"),
                                c("AOM","CONTROL"),
                                c("AOM","DSS"))


immuneMatrix =data.frame(GSE57533.immucc, Factor = samplesList$GROUP)
vars = colnames(immuneMatrix)[1:(dim(immuneMatrix)[2]-1)]
a = list()
testResults = data.frame()
for  (i in 1:length(vars)){
  
  data.formula = formula(paste(vars[i] ,"~ Factor"))
  testResults.temp =  compare_means(data.formula,comparisons = project.groupsToCompare, p.adjust.method = "fdr", method='wilcox.test',data = immuneMatrix)
  testResults = rbind(testResults,testResults.temp)
  
  #  pairwise.wilcox.test(x = immuneMatrix[,i], g =  immuneMatrix$Factor, p.adjust.method = "fdr")
  
  #esto muestra pvalue y no padjusted es una cagada
  p2 = ggboxplot(immuneMatrix,  x ="Factor",y = vars[i],
                 title = gsub("\\."," ",vars[i]),
                 subtitle = " ",
                 fill = "Factor", add = c("jitter"),width=0.5,
                 order = c( "CONTROL","DSS", "AOM"),
                 xlab = "",ylab ="TIL (%)",
                 repel = TRUE,combine = FALSE) +

    stat_compare_means(comparisons = project.groupsToCompare,
                       method = "wilcox.test",
                       p.adjust.methods = "fdr"
    ) 
  # +  stat_compare_means(label.y.npc = 1, label.x.npc = 0.3)  Kruskall
  
  #para agregar a mano padjusted, el tema es que te agrga todas las comapraciones y no las que queres  
  # stat_pvalue_manual(testResults, label = "p.adj",
  #                    y.position = seq(max(immuneMatrix[,i])+1, 10, length.out = dim(testResults)[1]))
  
  
  
  
  p2 = ggpar(p2, 
             font.title = c(14,"bold", "black"),palette = "npg"
  )
  p3 = p2
  p2 =  p2+rremove("x.text") +rremove("legend") +rremove("x.ticks")
  
  a[[i]] = p2
  #print(p2)  
  
  
}

HELPER_SAVE_DATA_FRAME(as.data.frame(testResults),"./results/tesis/ratón/immuCC/GSE57533/testResults.csv" )

p3  = p3 +rremove("x.title") +rremove("x.axis")+rremove("y.text") +rremove("y.ticks")+rremove("y.axis")+rremove("y.title") 
p3$labels$title = ""

p3 = ggpar(p3,
           legend = "right", legend.title = "Group", font.legend =c(12,"bold", "black"))
a[[i+1]] = p3

b =   plot_grid(plotlist =a,
                ncol = 5, nrow = 6)

ggsave(b,file=paste0("./results/tesis/ratón/immuCC/GSE57533/SVR_BoxplotWithPValues.pdf"), width = 15, height =25, limitsize = FALSE)
#*************************************************************************************************************************************
#*************************************************************************************************************************************


 #Humano
#1.  TCGA: primary tumor vs. solid tissue normal
# vs. 
# - GSE64658: AOM/DSS vs ctrl.
# - GSE57533: AOM/DSS vs ctrl.
# - GSE43338: AOM/DSS vs. ctrl.
# - GSE43338: Apc min/+ (spCRC) vs. ctrl.
heatmap.breakUp = 10
heatmap.breakDown = 7
heatmap.colorUp = colorRampPalette((brewer.pal(n = 8, name = "Reds")))(heatmap.breakUp)
heatmap.colorDown = colorRampPalette(rev(brewer.pal(n = 8, name = "Blues")))(heatmap.breakDown)
heatmap.colorZero = "white"
heatmap.colorNa = "white"
heatmap.colors=c( heatmap.colorDown,heatmap.colorZero,heatmap.colorUp)

humano = HELPER_LoadCSV("./results/tesis/humano/HumanoFinalMergeTesis.csv")
colnames(humano) = c("Genes",	
                    "Family",
                    "GSE43338: AOM/DSS VS Control",
                    "GSE43338: Apc min/+ (spCRC) VS Control",
                    "GSE57533: AOM/DSS VS Control",
                    "GSE57533: DSS VS Control",
                    "GSE64658: Late Distal AOM/DSS  VS Control",
                    "TCGA: Primary Tumor VS Solid Tissue Normal",
                    "TCGA: Stage I VS Solid Tissue Normal",
                    "TCGA: Stage II VS Solid Tissue Normal",
                    "TCGA: Stage III VS Solid Tissue Normal",
                    "TCGA: Stage IV VS Solid Tissue Normal",
                    "TCGA: Stage IV VS Stage I"	,
                    "TCGA: MSI-h VS Solid Tissue Normal")
humano = humano[,-6]
humano = humano[,1:7]
humano[,3:7][abs(humano[,3:7])<1.5] = NA
humano = na.omit(humano)
HELPER_SAVE_DATA_FRAME(humano,"./results/tesis/humano/humanoTumorVsCtrlVsRaton.csv")
rownames(humano) = humano$Genes
humano = humano[order(humano[,3],decreasing = TRUE),]
pheatmap(humano[,3:7], main = "TCGA  Vs Raton",
         cellheight = 10,cellwidth = 20,
         treeheight_row = 0, treeheight_col = 0, 
         gaps_col =  4,
         cluster_cols = FALSE,cluster_rows = FALSE,
         fontsize_row = 7,
         filename = "./results/tesis/humano/humanoTumorVsCtrlVsRaton.pdf")
#correlation
correlationResult = corr.test(x = humano[,3:7],use = "pairwise", method = "spearman", adjust = "BH", alpha=0.05, ci = FALSE)
HELPER_SAVE_DATA_FRAME(as.data.frame(correlationResult$r), "./results/tesis/humano/Correlation_R_humanoTumorVsCtrlVsRaton.csv")
HELPER_SAVE_DATA_FRAME(as.data.frame(correlationResult$p), "./results/tesis/humano/Correlation_PValue_humanoTumorVsCtrlVsRaton.csv")
pdf("./results/tesis/humano/CorrelationhumanoTumorVsCtrlVsStageVsRaton.pdf", width = 5, height= 5)
corrplot(corr = correlationResult$r,
         col = heatmap.colors,
         title ="Correlation", 
         diag = FALSE,
         type = "lower",
         tl.col="black", #Text label color and rotation
         number.cex =0.4,
         cl.cex = 0.35,
         tl.cex = 0.6,
         p.mat = correlationResult$p, sig.level = 0.05,insig = "blank",
         na.label = "  ",
         mar=c(0,0,2,0)
) 
dev.off()
#2
# TCGA: stage 1, 2, 3 y 4 tumor vs. solid tissue normal
# vs. 
# - GSE64658: AOM/DSS vs ctrl.
# - GSE57533: AOM/DSS vs ctrl.
# - GSE43338: AOM/DSS vs. ctrl.
# - GSE43338: Apc min/+ (spCRC) vs. ctrl


humano = HELPER_LoadCSV("./results/tesis/humano/HumanoFinalMergeTesis.csv")
colnames(humano) = c("Genes",	
                     "Family",
                     "GSE43338: AOM/DSS VS Control",
                     "GSE43338: Apc min/+ (spCRC) VS Control",
                     "GSE57533: AOM/DSS VS Control",
                     "GSE57533: DSS VS Control",
                     "GSE64658: Late Distal AOM/DSS  VS Control",
                     "TCGA: Primary Tumor VS Solid Tissue Normal",
                     "TCGA: Stage I VS Solid Tissue Normal",
                     "TCGA: Stage II VS Solid Tissue Normal",
                     "TCGA: Stage III VS Solid Tissue Normal",
                     "TCGA: Stage IV VS Solid Tissue Normal",
                     "TCGA: Stage IV VS Stage I"	,
                     "TCGA: MSI-h VS Solid Tissue Normal")
humano = humano[,-c(6,8,13,14)]

humano[,3:10][abs(humano[,3:10])<1.5] = NA
humano = na.omit(humano) #revisar por borra toda la linea si hay un solo NA
HELPER_SAVE_DATA_FRAME(humano,"./results/tesis/humano/humanoStageVsRaton.csv")
rownames(humano) = humano$Genes
humano = humano[order(humano[,3],decreasing = TRUE),]
pheatmap(humano[,3:10], main = "Tumor Stage Vs Raton",
         cellheight = 10,cellwidth = 20,
         gaps_col =  4,
         treeheight_row = 0, treeheight_col = 0, 
         cluster_cols = FALSE,cluster_rows = FALSE,
         fontsize_row = 7,
         filename = "./results/tesis/humano/humanoStageVsRaton.pdf")

#correlation
correlationResult = corr.test(x = humano[,3:10],use = "pairwise", method = "spearman", adjust = "BH", alpha=0.05, ci = FALSE)
HELPER_SAVE_DATA_FRAME(as.data.frame(correlationResult$r), "./results/tesis/humano/Correlation_R_humano_StageVsRaton.csv")
HELPER_SAVE_DATA_FRAME(as.data.frame(correlationResult$p), "./results/tesis/humano/Correlation_PValue_humano_StageVsRaton.csv")
pdf("./results/tesis/humano/Correlationhumano_StageVsRaton.pdf", width = 5, height= 8)
corrplot(corr = correlationResult$r,
         col = heatmap.colors,
         title ="Correlation", 
         diag = FALSE,
         type = "lower",
         tl.col="black", #Text label color and rotation
         number.cex =0.4,
         cl.cex = 0.35,
         tl.cex = 0.6,
         p.mat = correlationResult$p, sig.level = 0.05,insig = "blank",
         na.label = "  ",
         mar=c(0,0,2,0)
) 

dev.off()
#*************************************************************************************************************************************
#***OK REVISADO **********************************************************************************************************************
#*************************************************************************************************************************************
#3  Heatmap de DE de casi todos los experimentos
# Se lee una tabla solo con los que dieron siginificativo 
# Luego reemplazamos todos los que son Log2FC |<1.5| por NA 
# Elinamos los que tenga toda la fila con NA
# 
# - GSE64658: AOM/DSS vs ctrl.
# - GSE57533: AOM/DSS vs ctrl.
# - GSE43338: AOM/DSS vs. ctrl.
# - GSE43338: Apc min/+ (spCRC) vs. ctrl.
# - GSE57533: DSS vs. ctrl.
# - primary solid tumor vs. solid tissue normal
# - stage 1 vs. solid tissue normal
# - stage 2 vs. solid tissue normal
# - stage 3 vs. solid tissue normal
# - stage 4 vs. solid tissue normal
# - stage 4 vs. stage 1
# - MSI-h vs. solid tissue normal
#**************************************************************************************
allExperiments = HELPER_LoadCSV("./results/tesis/humano/HumanoFinalMergeTesis.csv")
colCount = dim(allExperiments)[2]
minFC2 = 1.5
removeMinFC = TRUE
fileName = "./results/Tesis/All.pdf"
colnames(allExperiments) = c("Genes",	
                     "Family",
                     "GSE43338: AOM/DSS VS Control",
                     "GSE43338: Apc min/+ (spCRC) VS Control",
                     "GSE57533: AOM/DSS VS Control",
                     "GSE57533: DSS VS Control",
                     "GSE64658: Late Distal AOM/DSS  VS Control",
                     "TCGA: Primary Tumor VS Solid Tissue Normal",
                     "TCGA: Stage I VS Solid Tissue Normal",
                     "TCGA: Stage II VS Solid Tissue Normal",
                     "TCGA: Stage III VS Solid Tissue Normal",
                     "TCGA: Stage IV VS Solid Tissue Normal",
                     "TCGA: Stage IV VS Stage I"	,
                     "TCGA: MSI-h VS Solid Tissue Normal")
title = paste0("Heatmap de Expresión Diferencial de genes de Glyco significativos.")

#Filtrado!!!!!!!! hay varias opciones
if (removeMinFC){
  #allExperiments[,3:colCount][abs(allExperiments[,3:colCount])<minFC2] = 0
  #allExperiments = allExperiments[rowSums(abs(allExperiments[3:14]),na.rm = TRUE )>0,] menos filtrado
  # allExperiments = allExperiments[rowSums(abs(allExperiments[3:7]),na.rm = TRUE )>0 & rowSums(abs(allExperiments[8:14]),na.rm = TRUE )>0,] #mas filtrado que este tanto como en ratón como en tcga
  #title = paste0("Heatmap de Expresión Diferencial de genes de Glyco significativos.\nSe removieron los Log2FC < |",minFC2,"| y luego las filas que no tenian valores")
  
  allExperiments = allExperiments[(rowSums(abs(allExperiments[3:7]),na.rm = TRUE )>0 & rowSums(abs(allExperiments[8:14]),na.rm = TRUE )==0) | (rowSums(abs(allExperiments[3:7]),na.rm = TRUE )==0 & rowSums(abs(allExperiments[8:14]),na.rm = TRUE )>0),] #Que este Solo en TCGA o Solo en Raton
  title= paste0("Heatmap de Expresión Diferencial de genes de Glyco significativos.\nSe removieron los Log2FC < |",minFC2,"|. Se expresa solo en Ratón o solo en Humano")
  
}
allExperiments$Genes = droplevels(allExperiments$Genes)
allExperiments$Family = droplevels(allExperiments$Family)
HELPER_SAVE_DATA_FRAME(allExperiments,"./results/Tesis/All.csv")

#escala
heatmap.breaks = unique(c(
  seq(min(allExperiments[,3:14], na.rm = TRUE),-0.5, length=heatmap.breakDown),
  0,
  rev(seq(max(allExperiments[,3:14], na.rm = TRUE),0.5, length=heatmap.breakUp))
))
legendBreaks = floor(heatmap.breaks)


#rowGaps = as.data.frame(table(allExperiments$Family))


 dim( allExperiments)
#allExperiments = allExperiments[rowSums(abs(allExperiments[3:14]))>0,]
 
 annotationRow =data.frame(Family = allExperiments$Family)
 rownames(annotationRow) = allExperiments$Genes
 annotationCol =data.frame(Type = c("Mouse","Mouse","Mouse","Mouse","Mouse", 
                                    "TCGA","TCGA","TCGA","TCGA","TCGA","TCGA","TCGA"))
 rownames(annotationCol) = colnames(allExperiments)[3:14]
 
 type.colors = settings.graphics.colors[1:2]
 names(type.colors) = c("Mouse", "TCGA")
 
 family.colors = settings.graphics.colors[3:(length(levels(allExperiments$Family))+2)]
 names(family.colors) = levels(allExperiments$Family)
 annotationColor = list(Type = type.colors,
                        Family = family.colors)

 #separación por familia de enzimas
rowGaps.freq = as.data.frame(table(allExperiments$Family))
rowGaps.freq = rowGaps.freq$Freq
rowGaps = 0
for (i in 1:length(rowGaps.freq)){
  rowGaps = rowGaps+ rowGaps.freq[i]
  rowGaps.freq[i] = rowGaps
}
rowGaps.freq = rowGaps.freq[1:(length(rowGaps.freq)-1)]
rownames(allExperiments) = allExperiments$Genes
allExperiments = allExperiments[order(allExperiments$Family,allExperiments$Genes),]
pheatmap(allExperiments[,3:14], main = title,
         cellheight = 10,cellwidth = 30,
         treeheight_row = 120, treeheight_col = 30, 
         color=heatmap.colors,breaks = heatmap.breaks,
         gaps_col =  5,
         gaps_row = rowGaps.freq,
         legend_breaks = legendBreaks,
         display_numbers = FALSE,
         annotation_colors = annotationColor,
         annotation_row = annotationRow,
         annotation_col = annotationCol,
         cluster_cols = FALSE,cluster_rows = FALSE,
         fontsize_row = 7,
         filename = fileName)
#*************************************************************************************************************************************
#*************************************************************************************************************************************
#*************************************************************************************************************************************

#correlation
correlationResult = corr.test(x = humano[,3:14],use = "pairwise", method = "spearman", adjust = "BH", alpha=0.05, ci = FALSE)
HELPER_SAVE_DATA_FRAME(as.data.frame(correlationResult$r), "./results/tesis/humano/Correlation_R_All.csv")
HELPER_SAVE_DATA_FRAME(as.data.frame(correlationResult$p), "./results/tesis/humano/Correlation_PValue_All.csv")
pdf("./results/tesis/humano/Correlation_All.pdf", width = 5, height= 8)
corrplot(corr = correlationResult$r,
         col = heatmap.colors,
         title ="Correlation", 
         diag = FALSE,
         type = "lower",
         tl.col="black", #Text label color and rotation
         number.cex =0.4,
         cl.cex = 0.35,
         tl.cex = 0.6,
         p.mat = correlationResult$p, sig.level = 0.05,insig = "blank",
         na.label = "  ",
         mar=c(0,0,2,0)
) 

dev.off()
#*************************************************************************************************************************************
#*************************************************************************************************************************************
#*************************************************************************************************************************************
#4. Heatmap de expresion solo de genes de glyco

#*************************************************************************************************************************************
#GSE57533
#*************************************************************************************************************************************

geneList  = HELPER_LoadCSV("Resources/GlycoGeneList/Glyco_MM_GeneList.csv")
samples = HELPER_LoadTSV("Results/MOUSE/GSE57533/Samples/GSE57533_SampleList.csv")
samples$SAMPLE = str_replace(samples$SAMPLE,"_"," ")
GSE57533.countMatrix = HELPER_LoadTSV("./Results/MOUSE/GSE57533/HTSEQ/CountMatrix/rawCountMatrixAnnotated.csv")
colnames(GSE57533.countMatrix)[2:11] = samples$SAMPLE

#PCA
GSE57533.pca = DA_GetPCA(t(GSE57533.countMatrix[,2:11]))
fviz_pca_ind(GSE57533.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

DA_PlotPCA(GSE57533.pca)


#Distance Maxtrix
DA_PlotSamplesDistance(as.matrix(t(GSE57533.countMatrix[,2:11])),sampleNames = samples$SAMPLE, sampleGroups= samples$GROUP,size = 20, filePath = "./results/Tesis/GSE57533_SampleDistance.pdf")


#Glycon Gene filter
GSE57533.countMatrix= merge(GSE57533.countMatrix,geneList, by.x = "ENSENMBL_GENE_ID", by.y = "ensembl_gene_id")

#escalamos para normalizar con log2
#GSE57533.countMatrix[,2:11]= log2(GSE57533.countMatrix[,2:11]+0.1)



annotationRow =data.frame(Family = GSE57533.countMatrix$Group)
rownames(annotationRow) = GSE57533.countMatrix$external_gene_name.x
annotationCol =data.frame(Type = samples$GROUP)
rownames(annotationCol) = samples$SAMPLE

type.colors = settings.graphics.colors[1:length(levels(factor(samples$GROUP)))]
names(type.colors) = levels(factor(samples$GROUP))

family.colors = settings.graphics.colors[3:(length(levels(GSE57533.countMatrix$Group))+2)]
names(family.colors) = levels(GSE57533.countMatrix$Group)
annotationColor = list(Type = type.colors,
                       Family = family.colors)


#GSE57533.countMatrix = GSE57533.countMatrix[order(GSE57533.countMatrix$Group,GSE57533.countMatrix$external_gene_name.x),]

fileName = "./results/Tesis/GSE57533_GlycoCountMatrixHeatmap.pdf"
title = "GSE57533   Glyco Genes Count Matrix Expression \nZ-Score normalization"
pheatmap(GSE57533.countMatrix[,2:11], main = title,
         scale = "row",
        # cellheight = 1,
        cellwidth = 30,
         treeheight_row = 120, treeheight_col = 20, 
         color=heatmap.colorUp,#breaks = heatmap.breaks,
       #  gaps_col =  rowGaps.col,
      #   gaps_row = rowGaps.freq,
        # legend_breaks = legendBreaks,
         display_numbers = FALSE,
         annotation_colors = annotationColor,
        # annotation_row = annotationRow,
         annotation_col = annotationCol,
      show_rownames = FALSE,
         cluster_cols = TRUE,cluster_rows = FALSE,
         #fontsize_row = 7,
         filename = fileName)

pdf("./results/Tesis/GSE57533_GlycoCountMatrixHeatmap2.pdf")
heatmap(as.matrix(GSE57533.countMatrix[,2:11]))
dev.off()
#*************************************************************************************************************************************
#*************************************************************************************************************************************

#esto viejo no va
DESEQ2_volcano(DESeqResult,
               cutoff = alpha,
               title = title ,
               minLog2FC = log2FC, 
               colors = c("gray", "yellow",rgb(102/255,151/255,234/255,1),rgb(176/255,36/255,40/255,1) ),
               xMaxValues = c(x.min,x.max),
               yMaxValues = c(0,y.max),
               showLegend = FALSE,
               geneLabelMinLog2FC = glycoGenes.log2FC,
               geneLabelMinAlpha = glycoGenes.alpha,
               showGeneNames = TRUE,
               showGeneNameList = glycoGenes.toShow,
               showlines = TRUE,
               labelOffset = 0.25,
               drawLabelBackground = TRUE,
               filePath =  fileName)


#*************************************************************************************************************************************
#Volcano plot
#*************************************************************************************************************************************
#GSE57533
#*************************************************************************************************************************************
source("./src/shared/functions/DEA_Function.r")
source("./src/shared/functions/HELPER_Function.r")
library(EnhancedVolcano)
geneList  = HELPER_LoadCSV("Resources/GlycoGeneList/Glyco_MM_GeneList.csv")
#geneList  = HELPER_LoadCSV("Resources/GlycoGeneList/Glyco_HS_GeneList.csv")

fileName = "./../Volcanoplot"
title = "colorectal_control_epithelium_colitis_associated_VS_colorectal_control_epithelium_sporadic"
data = ""
data = HELPER_LoadTSV("./Results/MOUSE/GSE43338/DEA/colorectal_control_epithelium_colitis_associated_VS_colorectal_control_epithelium_sporadic/DEA/RES_ALL_GENES.csv")
View(data)
#data = HELPER_LoadTSV("./Results/HUMAN/GDC/TCGA/COAD_TUMOR_MSI/HTSEQ/DEA/mss_VS_Solid Tissue Normal/DEA/RES_ALL_GENES.csv")


data = data[data$padj != 1,]
log2FC = 1
alpha = 0.05
glycoGenes.alpha = 0.05
glycoGenes.log2FC = 2
table (data$padj<alpha)
table (data$padj<alpha  & tolower(data$external_gene_name) %in% tolower(geneList$external_gene_name))
table (data$padj<alpha & abs(data$log2FoldChange) >log2FC)
table (data$padj<alpha & abs(data$log2FoldChange) >log2FC & tolower(data$external_gene_name) %in% tolower(geneList$external_gene_name))
table (data$padj<alpha & abs(data$log2FoldChange) >glycoGenes.log2FC)
table (data$padj<alpha & abs(data$log2FoldChange) >glycoGenes.log2FC & tolower(data$external_gene_name) %in% tolower(geneList$external_gene_name))

glycoGenes.toShow = data$external_gene_name[data$padj<alpha & abs(data$log2FoldChange) >glycoGenes.log2FC & tolower(data$external_gene_name) %in% tolower(geneList$external_gene_name)]
length(glycoGenes.toShow)
#remove actb
glycoGenes.toShow = glycoGenes.toShow[-3] 

HELPER_SAVE_DATA_FRAME(data.frame(glycoGenes.toShow),"./../GSE57533_CCRAC_VS_Control.tsv")
length(glycoGenes.toShow)
x.min =  min(data$log2FoldChange,na.rm = TRUE) -1
x.max = max(data$log2FoldChange,na.rm = TRUE) 
y.max = 7
DESeqResult = data[!is.na(data$padj),] 

keyvals.colour <- 
  ifelse(DESeqResult$padj > alpha, "gray", 
         ifelse(DESeqResult$padj < alpha &DESeqResult$log2FoldChange > -log2FC & DESeqResult$log2FoldChange < log2FC, "yellow",
                ifelse(DESeqResult$padj < alpha &DESeqResult$log2FoldChange < -log2FC, rgb(102/255,151/255,234/255,1),
                       ifelse(DESeqResult$padj < alpha & DESeqResult$log2FoldChange > log2FC, rgb(176/255,36/255,40/255,1),
                              'gray'))))

keyvals.colour[is.na(keyvals.colour)] <- 'gray'
names(keyvals.colour)[keyvals.colour == 'gray'] <- 'No significativo'
names(keyvals.colour)[keyvals.colour == rgb(176/255,36/255,40/255,1)] <- 'Up'
names(keyvals.colour)[keyvals.colour == 'yellow'] <- '|Log2FC| <1 '
names(keyvals.colour)[keyvals.colour == rgb(102/255,151/255,234/255,1)] <- 'Down'



pdf(paste0(fileName,".pdf"), width=3.5, height=6)
EnhancedVolcano(DESeqResult,
                lab = DESeqResult$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = glycoGenes.toShow,
                title = title,
                titleLabSize = 4,
                subtitle = "",
                captionLabSize = 0,
                xlim = c(x.min,x.max),
                ylim = c(0,y.max),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~'(p-adj)'),
                pCutoff = alpha,
                FCcutoff = log2FC,
               # col = c("gray","gray", "yellow",rgb(176/255,36/255,40/255,1) ),
               colCustom = keyvals.colour,
                boxedlabels = TRUE,
                colAlpha = 0.4,
               axisLabSize = 12,
                legendPosition = 'bottom',
                legendLabSize = 0,
                legendIconSize = 0,
               
               legendVisible = FALSE,
                drawConnectors = TRUE,

                widthConnectors = 0.3,
                colConnectors = 'black')
dev.off()
#**********************
# Media de infitlrado en xcell
#**********************
library(readr)
library(psych)
XCELL <- read_delim("/Volumes/Externo/Google Drive/Bioinformática/PFI/Tesina/Tesis Final/Resultados Finales/Human/GDC/TCGA/COAD/COAD_TUMOR_STAGE/immune/xCell/XCELL.csv", 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

COAD_TUMOR_STAGE_SampleList <- read_delim("Results/Human/GDC/TCGA/COAD_TUMOR_STAGE/Samples/COAD_TUMOR_STAGE_SampleList.csv", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE)
View(COAD_TUMOR_STAGE_SampleList)
View(XCELL)
tablaPromedios = data.frame(Grupo= COAD_TUMOR_STAGE_SampleList$GROUP, Estadio = COAD_TUMOR_STAGE_SampleList$tumor_stage,MSI = COAD_TUMOR_STAGE_SampleList$msi_status, dMMR = COAD_TUMOR_STAGE_SampleList$loss_expression_of_mismatch_repair_proteins_by_ihc,XCELL[,2:(dim(XCELL)[2]-3)] )
View(tablaPromedios)

library(dplyr)
res = tablaPromedios %>% 
  group_by(Grupo) %>%
  summarise(across(everything(), mean))

res = t(res)
View(res)
View(describeBy(tablaPromedios[5:dim(tablaPromedios)[2]], tablaPromedios$Grupo, mat = TRUE) )
