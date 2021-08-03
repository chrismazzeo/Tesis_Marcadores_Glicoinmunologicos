#TODO !!!! hay que sacar duntest y hacer wilcoxtest 
#pairwise.wilcox.test(x = immuneMatrix[,1], g =  immuneMatrix$Factor, p.adjust.method = "fdr")


IMMUNE_ImmuCC_RNASEQ = function (immuCCPath = "./ExternalLibs/ImmuCC/RNASeq/",rawExpressionMatrix,outputPath){
  source(paste0(immuCCPath,"MouseHTSeq_counts_stat.r"))

  immuCC = getRNASEQInfiltrations(rawExpressionMatrix,paste0(immuCCPath,"receptor.ensemble.merge.RData"))

  immuCC =  HELPER_RowNamesAsFirstColumn(immuCC, "ENSENMBL_GENE_ID")
  HELPER_SAVE_DATA_FRAME(immuCC,outputPath)
  return (immuCC)
}


IMMUNE_ImmuCC_RNASEQ_STATS = function (filePath, outputDir, title,samplesName,group,decovolutionalMethod ="SVR",groupToCompare){
  
  dir.create(file.path( outputDir), recursive = TRUE, showWarnings = FALSE)
  
  immuCCResult = HELPER_LoadCSV(filePath)
  immuCCResult = immuCCResult[,-1]
  rownames(immuCCResult) = samplesName
  
  
  immuCCResult.save = HELPER_RowNamesAsFirstColumn(immuCCResult, "SAMPLE")
  HELPER_SAVE_DATA_FRAME(immuCCResult.save,paste0(outputDir,decovolutionalMethod,".csv"))
  

  rowAnnotation <- data.frame(GROUP = group)
  rownames(rowAnnotation) = rownames(immuCCResult)
  
  pheatmap(as.matrix(immuCCResult),display_numbers = FALSE, main = title,
            cluster_rows = TRUE, cluster_cols = FALSE, 
            cellheight = 10,cellwidth = 20,
            treeheight_row = 0, treeheight_col = 0,
            annotation_row = rowAnnotation , 
            fontsize_row = 7,
            filename = paste0(outputDir,decovolutionalMethod,"_heatmap.pdf"))

  if (length(levels(group))>1){
    IMMUNE_BOX_PLOT_WITH_PVALUES(immuneMatrix = immuCCResult,
                                 groups = group,
                                 outputDir = outputDir,
                                 title ="",
                                 decovolutionalMethod = decovolutionalMethod,
                                 width = 20+  length(groupToCompare)*0.1, 
                                 height = 0.55*dim(immuneMatrix)[2]+5+ length(groupToCompare)*0.5,
                                 groupToCompare = groupToCompare)
    IMMUNE_CELL_COMPARE(data = immuCCResult,
                        group = group,
                        outputDir= outputDir,
                        decovolutionalMethod = decovolutionalMethod
    )
  }
  
  
  
  IMMUNE_STACKED_BAR_PLOT(
    immuneMatrix = immuCCResult,
    groups = group,
    outputDir = outputDir,
    title = "All Samples",
    xTitle ="", 
    yTitle = "(%) Immune Infiltration Cell Type" ,
    leftOffset = 10,
    graphicsWidth = 3,
    graphicsHeight = 5,
    decovolutionalMethod = decovolutionalMethod
  )
               
 
  IMMUNE_CELL_CORRELATION(immuneMatrix = immuCCResult,
                          group = group,
                          outputDir = outputDir,
                          decovolutionalMethod = decovolutionalMethod, 
                          padj = 0.05,
                          method = "spearman")
  
  IMMUNE_CELL_COMPARE(data = immuCCResult,
                      group = group,
                      outputDir= outputDir,
                      decovolutionalMethod = decovolutionalMethod
  )
  
  #*******************************************************************************************************
  #DA
  #*******************************************************************************************************


  immuCC.tsne.2d2 = DA_GetTSNE(dataFrame =immuCCResult, dimensions = 2,perplexity = 2, checkDuplicates = FALSE)
  DA_Plot2D(immuCC.tsne.2d2$Y, title = "TSNE 2D ImmuCC", colorsArray = settings.color, group = group,outputDir = outputDir, fileName = paste0(decovolutionalMethod,"_TSNE2d.pdf"))
  
  #immuCC.PCA.svr = DA_GetPCA(immuCCResult, scale = TRUE)
 # DA_PlotPCA(immuCC.PCA.svr)
  
  return (immuCCResult)
}


#***************************************************************
#data frame rowsnames = samples, colnames = only cell type, is the results of Cibersort or TIMER | remove other columns 
#group = factor

#***************************************************************

IMMUNE_STACKED_BAR_PLOT = function (immuneMatrix,groups, outputDir,title="", xTitle ="", yTitle ="" , colorArray, groupColorArray, leftOffset = 10, graphicsWidth = NA, graphicsHeight = NA,decovolutionalMethod = "SVR"){
 immuneMatrix.factor =data.frame(immuneMatrix, Factor = groups)

   #grouped -stacked barplot
  d <- melt(immuneMatrix.factor, id.vars="Factor")
  b =    aggregate(d$value, list(d$Factor,d$variable), mean)
  colnames (b) = c("Factor","Type", "Value")
 
  p4 = ggplot(data=b, aes(x=Factor, y=Value, fill = Type))  + geom_bar(stat="identity",width = 0.5)+
    ggtitle("Grouped Barplot")+ 
    xlab("")+
    ylab(" (%) Immune Infiltration Cell Type") +
    settings.graphics.fillColor+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text = element_text(size=7),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
 
  p4
  ggsave(p4,file=paste0(outputDir,decovolutionalMethod,"_Grouped-StackedBarplot.pdf"), width =  length(levels(groups))*0.5 + 0.15*dim(immuneMatrix)[1],
         height = dim(immuneMatrix)[1]*.08, limitsize = FALSE)

  
  #only grouped  
  immuneMatrix =t (immuneMatrix)
  print (head(immuneMatrix))
  d <- melt(immuneMatrix, id.vars="Factor")
  colnames (d) = c("Type", "Sample", "Value")
  p5<-ggplot(data=d, aes(x=Sample, y=Value, fill = Type)) +
    geom_bar(stat="identity") +
  
    ggtitle("Grouped Barplot")+ 
    ylab(" (%) Immune Infiltration Cell Type") +
    settings.graphics.fillColor+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
  
  p5
  ggsave(p5,file=paste0(outputDir,decovolutionalMethod,"_StackedBarplot.pdf"), width = dim(immuneMatrix)[2]*0.3 +0.15*dim(immuneMatrix)[1], height =  dim(immuneMatrix)[1]*.23, limitsize = FALSE)
  
}
GET_COMPARE_LIST = function(groups){

  compareList = list()
  pos = 1
  items = levels(groups)
  for (i in 1:(length(items)-1)){
    rest = items[(i+1):(length(items))]
    for (x in 1:length(rest)){
      compareList[[pos]] = c(as.character(items[i]), as.character(rest[x]))
      pos = pos+1
    }
  }
 return (compareList)
  
}
IMMUNE_BOX_PLOT_WITH_PVALUES =function (immuneMatrix,groups,outputDir, title ="",ylab =" (%) Immune Infiltration Cell Type",decovolutionalMethod = "SVR", width =NA, height = NA, groupToCompare = NA){
  immuneMatrix =data.frame(immuneMatrix, Factor = groups)
  print (head(immuneMatrix))

  #boxplot
  p1 = ggboxplot(immuneMatrix,  x ="Factor",y = colnames(immuneMatrix)[1:(dim(immuneMatrix)[2]-1)],merge = TRUE, ylab = ylab )+
    settings.graphics.theme +
   
    settings.graphics.borderColor
  p1
  p1 =  ggpar(p1, legend = "right") 
  ggsave(p1,file=paste0(outputDir,decovolutionalMethod,"_Boxplot.pdf"), width  = 15, 
         height =4, limitsize = FALSE)
  
 
#boxplot with pvalues
  height.dist = 20
  height.Offset = 110
  if (is.na(groupToCompare)){
    my_comparisons <- GET_COMPARE_LIST(groups)
  }
  else{
    my_comparisons = groupToCompare
  }
  
  
  p2 = ggboxplot(immuneMatrix,  x ="Factor",y = colnames(immuneMatrix)[1:(dim(immuneMatrix)[2]-1)],combine = TRUE,repel = FALSE,
                 color = "Factor", add = c("boxplot"),width=1,
                 xlab = TRUE,ylab =ylab)+
    stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y =seq(height.Offset,height.Offset+length(my_comparisons)*height.dist,height.dist))+
    stat_compare_means(label.y = height.Offset+20 + length(my_comparisons) *height.dist)+
    settings.graphics.theme +
    settings.graphics.borderColor
  p2
  ggsave(p2,file=paste0(outputDir,decovolutionalMethod,"_BoxplotWithPValues.pdf"), width = width, height = height*1.5, limitsize = FALSE)
  
  
  # both
  p3 = grid.arrange(p2, p1, nrow = 2)
  p3
  ggsave(p3,file=paste0(outputDir,decovolutionalMethod,"_BoxplotWithPValues2.pdf"), width = width, height = height, limitsize = FALSE)

}

IMMUNE_CELL_CORRELATION = function (immuneMatrix,groups,outputDir, decovolutionalMethod = "SVR",padj, method = "spearman"){
  #***************************************************************
  #correlacion x grupos


  #***************************************************************

  for (i in 1: length(levels(factor(groups)))){
    corMatrix = immuneMatrix[HELPER_ConvertFactorToInteger(groups) == i,]
    #si el the standard deviation is zero no saco nada
    
      correlationResult = corr.test(x = corMatrix,use = "pairwise", method = method, adjust = "BH", alpha=padj, ci = FALSE)
    
      cor = HELPER_RowNamesAsFirstColumn(correlationResult$r, "CEL")
      padj = HELPER_RowNamesAsFirstColumn(correlationResult$p, "CEL")
      HELPER_SAVE_DATA_FRAME(cor,paste0(outputDir,decovolutionalMethod,"_Correlation_Group_",paste0(levels(factor(groups))[i],".csv")))
      HELPER_SAVE_DATA_FRAME(padj,paste0(outputDir,decovolutionalMethod,"_Correlation_PAdjvalue_Group_",paste0(levels(factor(groups))[i],".csv")))

      pdf(file=paste0(outputDir,decovolutionalMethod,"_Correlation_Group_",paste0(levels(factor(groups))[i],".pdf")))
  
      corrplot(corr = correlationResult$r,
             title = paste0("Correlation ",levels(factor(groups))[i]), 
             diag = FALSE,
             type = "lower",
             tl.col="black", #Text label color and rotation
             number.cex =0.7,  addCoef.col = "white",addgrid.col   = "gray",
             p.mat = correlationResult$p, sig.level = padj, insig = "label_sig",
             na.label = "  ",
             mar=c(0,0,1,0)
      ) 
      dev.off()
      corrplot(corr = correlationResult$r,
             title = paste0("Correlation ",levels(factor(groups))[i]), 
             diag = FALSE,
             type = "lower",
             tl.col="black", #Text label color and rotation
             number.cex =0.7,  addCoef.col = "white",addgrid.col   = "gray",
             p.mat = correlationResult$p, sig.level = padj, insig = "blank",
             na.label = "  ",
             mar=c(0,0,1,0)
      )
    
  }
}
IMMUNE_CELL_COMPARE = function(data, group,decovolutionalMethod,  outputDir, adjMethod = "bonferroni"){

  numbersOfImmneCells = ncol(data)
  data = data.frame(data, Factor = factor(group))
  
  #***************************************************************
  #vemos normalidad multivariada x grupo 
  #***************************************************************
  normality = TRUE
  sink(paste0(outputDir,decovolutionalMethod, "_mvShapiroResults.txt"))
  for (i in 1:length(levels(data$Factor))){
      tempShapiroResult <- try(mvShapiro.Test(as.matrix(data[data$Factor == levels(data$Factor)[i],1:numbersOfImmneCells])))
      if(inherits(tempShapiroResult, "try-error"))
      {
        #error handling code, maybe just skip this iteration using
        print (paste0("Group ",levels(data$Factor)[i], " Not enougth samples to run test"))
        normality = FALSE
        next
      }
      else{
        if (tempShapiroResult$p.value >0.05){
          print (paste0("Group : ",levels(data$Factor)[i], "Normality TRUE"))
        }
        else{
          print (paste0("Group : ",levels(data$Factor)[i], "Normality FALSE"))
          normality = FALSE
        }
        print (tempShapiroResult)
      }
  }
  sink(type = "message")
  sink()
  if (!normality){
          #***************************************************************
          #comparacion entre los  grupos NO Parametrico 
          #Para cada variablie kruskall --> si da significativo comapracion por grupos y correccion por multiple corrcciones 
          #dunntest hace kruskall + comapracion de a 2 con correccion
          #***************************************************************
          
    #combinatoria para saber la cantidad de columnas
          cols =  factorial(length(levels(data$Factor)))/(2*(factorial (length(levels(data$Factor)) - 2)))
          cols = cols +1 #para guardar el resultado de Kruskall
                   
          #preparamos el dataframe para guardar la lista
          testGroupMatrix.results = data.frame( matrix(nrow = numbersOfImmneCells,ncol =cols))
          rownames(testGroupMatrix.results) =colnames(data)[1:numbersOfImmneCells]
          colnames(testGroupMatrix.results)[1] = "Kruskall-Wallis"
          
          
          #comparamos
          for (i in 1:(numbersOfImmneCells)){
            #hacemos kruskall para ver si hay diferencias entre grupos
            kuskallResult = kruskal.test(data[,i] ,g  = data$Factor)
            testGroupMatrix.results[i,1]= kuskallResult$p.value
            
            #hacemos comparacion entre grupos si kruskall da significativo
            if (!is.na(kuskallResult$p.value) & kuskallResult$p.value <0.05){
                tempResults = try(dunnTest(data[,i] ,g  = data$Factor,method="bonferroni", two.sided = TRUE))
                if(inherits(tempResults, "try-error"))
                {
                  #error handling code, maybe just skip this iteration using
                  
                }
                else{
                  for (j in 1:dim(tempResults$res)[1]){
                    testGroupMatrix.results[i,1+j] = tempResults$res$P.adj[j]
                    colnames(testGroupMatrix.results)[1+j] = as.character(tempResults$res$Comparison[j])
                  }
                }
                
            }
          }
  }
  testGroupMatrix.results
  
  
  #***************************************************************
  # Tabla mean de cada variable por grupo + adjusted p value
  #TODO guardar la sumatoria de cada columna, debe dar 100%, agregar "%" al titulo, agregar columna de desvio
  #***************************************************************
  resulTable = data.frame(Cell = rownames(testGroupMatrix.results), padj = testGroupMatrix.results)
  write_csv(resulTable, path = paste0(outputDir,decovolutionalMethod,"_Resultados_Comapración_entre_grupos.csv"))
}

#***************************************************************
# give the correct format to immune matrix
# type= 0.Cibersort  1. Timer
#***************************************************************
replaceRowsNames = function (data, column){
  data = as.data.frame(data)
  rownames(data) = data[,column]
  
  data =data[,-column]  
  return (data)
}


#*******
# testGroupMatrix.results %>%
#   mutate(
#     Var = row.names(.)
#   )%>%
#   mutate_if(is.numeric, function(x) {
#     x= ifelse(x < pValueCutOff,
#               cell_spec(x, color = "red", bold = T),
#               cell_spec(x, color = "black", italic = T))
#   })%>%
#   kable(escape = F, align = "c",caption  = "Adjusted p Value") %>%
#   kable_styling(bootstrap_options = c("striped", "hover", "condensed", full_width = T))  %>%
#   cat(., file = output.infiltration.cibersort.groupTest)




IMMUNE_Humman_RNASEQ_STATS = function (immuneMatrix,pvalue, outputDir, title,group,decovolutionalMethod ="CIBERSORT", groupToCompare, orderFactor = NA, width = 15, height= 20, colCount ){
 # immuCCResult = HELPER_LoadCSV(filePath)
 # immuCCResult = immuCCResult[,-1]
#  rownames(immuCCResult) = samplesName
  dir.create(file.path( outputDir), recursive = TRUE, showWarnings = FALSE)
  
  immuCCResult = as.data.frame(immuneMatrix)
  
  immuCCResult.save = HELPER_RowNamesAsFirstColumn(immuCCResult, "SAMPLE")
  HELPER_SAVE_DATA_FRAME(immuCCResult.save,paste0(outputDir,decovolutionalMethod,".csv"))
  
  if(decovolutionalMethod == "CIBERSORT"){
    group = group[immuCCResult$`P-value`<pvalue]
    group = droplevels(group)
    samplesList.inmune.filters = samplesList[immuCCResult$`P-value`<pvalue,]
    immuCCResult = immuCCResult[immuCCResult$`P-value`<pvalue,1:22]
    
  }
  else{
    samplesList.inmune.filters = samplesList
  }
  
  if(decovolutionalMethod == "XCELL"){
    immuCCResult = immuCCResult[,1:64]
    
  }
  if(decovolutionalMethod != "SVR"){
    immuCCResult = immuCCResult*100
    
  }
  
  
  rowAnnotation <- data.frame(GROUP = group)
  rownames(rowAnnotation) = rownames(immuCCResult)
  
  pheatmap(as.matrix(immuCCResult), main = title, cluster_rows = TRUE, cluster_cols = FALSE, 
           annotation_row = rowAnnotation ,  treeheight_row = 0,
           width =(0.37*dim(immuCCResult)[2]+5) ,height = (0.15*dim(immuCCResult)[1]+5),
           filename = paste0(outputDir,decovolutionalMethod,"_heatmap.pdf"))
  
  
  #heatmap in reds
  
  rowAnnotation2 <- data.frame(Type = samplesList.inmune.filters$GROUP, Stage = samplesList.inmune.filters$tumor_stage, MSI = samplesList.inmune.filters$msi_status, Lost_Repair = samplesList.inmune.filters$loss_expression_of_mismatch_repair_proteins_by_ihc)
  rownames(rowAnnotation2) = rownames(immuCCResult)
  rowAnnotation2$Stage[rowAnnotation2$Type == "Solid Tissue Normal"] = NA
  rowAnnotation2$Stage =droplevels(rowAnnotation2$Stage)
  
  rowAnnotation2$MSI [rowAnnotation2$Type == "Solid Tissue Normal"] = NA
  rowAnnotation2$MSI =droplevels(rowAnnotation2$MSI)
  
  rowAnnotation2$Lost_Repair [rowAnnotation2$Type == "Solid Tissue Normal"] = NA
  rowAnnotation2$Lost_Repair =droplevels(rowAnnotation2$Lost_Repair)
  
  levels(rowAnnotation2$Type) = c("Tumor", "Control")
  levels(rowAnnotation2$Stage)   = gsub(" ","",levels(rowAnnotation2$Stage)  )
  levels(rowAnnotation2$MSI)  = gsub("-","_",levels(rowAnnotation2$MSI)  )
  levels(rowAnnotation2$Lost_Repair)
  ann_colors = list(
     Type = c(Tumor = "#4f50ff", Control = "#5a645f"),
     Stage =  c(stagei = "#749b57", stageii = "#efe685", stageiii = "#466983", stageiv = "#5db0dd","white" ),
     MSI = c(msi_h ="#802168", msi_l = "#6cd66b" ,mss = "#d494a7","white"),
     Lost_Repair = c(yes = "#924720", no = "#837a8d","white")
  )
  
  colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
 # pheatmap(as.matrix(immuCCResult), main = title, cluster_rows = TRUE, cluster_cols = TRUE, 
  #         annotation_row = rowAnnotation2 ,  treeheight_col = 0,
   #        col=colors,annotation_colors = ann_colors,
    #       width =(0.37*dim(immuCCResult)[2]+5) ,height = (0.15*dim(immuCCResult)[1]+5),
     #      filename = paste0(outputDir,decovolutionalMethod,"_heatmap2.pdf"))
  
# IMMUNE_STACKED_BAR_PLOT(
#  immuneMatrix = immuCCResult,
#    groups = group,
#    outputDir = outputDir,
#    title = "All Samples",
#    xTitle ="", 
#    yTitle = "(%) Immune Infiltration Cell Type" ,
#    leftOffset = 10,
#    graphicsWidth = 40,
#    graphicsHeight = 12,
#    decovolutionalMethod = decovolutionalMethod
#  )
  
  
 # IMMUNE_CELL_CORRELATION(immuneMatrix = immuCCResult,
  #                        group = group,
  #                        outputDir = outputDir,
  #                       decovolutionalMethod = decovolutionalMethod, 
  #                       padj = 0.05,
  #                    method = "spearman")
  if (length(levels(group))>1){
    BoxPlotWithPvalues(immuneMatrix = immuCCResult,
                                 groups = group,
                                 outputDir = outputDir,
                                 title ="",
                                 decovolutionalMethod = decovolutionalMethod,
                                 width = width, 
                                 height = height,
                                 orderFactor = orderFactor,
                                 groupToCompare = groupToCompare,
                                colCount = colCount)
    IMMUNE_CELL_COMPARE(data = immuCCResult,
                        group = group,
                        outputDir= outputDir,
                        decovolutionalMethod = decovolutionalMethod
    )
  }
  #*******************************************************************************************************
  #DA
  #*******************************************************************************************************
  
#  try({
#    immuCC.tsne.3d2 = DA_GetTSNE(dataFrame =immuCCResult, dimensions = 3,perplexity = 2, checkDuplicates = FALSE)
#    DA_Plot3D(immuCC.tsne.3d2$Y, title = "TSNE 3D ImmuCC", colorsArray = settings.color, group = group,outputDir = outputDir, fileName = paste0(decovolutionalMethod,"_TSNE3d.pdf"))
#  })
 
   #immuCC.PCA.svr = DA_GetPCA(immuCCResult, scale = TRUE)
  # DA_PlotPCA(immuCC.PCA.svr)
  
#  res.pca = princomp(immuCCResult, cor = FALSE, scores = TRUE)
#  pdf(paste0(outputDir,"pca.pdf"))
#  print(fviz_eig(res.pca))
#  print(fviz_pca(res.pca,     # Avoid text overlapping,
#           label = "var",
#           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#           col.var = "contrib",
#           repel = TRUE ))
 # dev.off()
  
  
  return (immuCCResult)
}

BoxPlotWithPvalues = function(immuneMatrix,groups,outputDir, title ="",ylab =" (%) Immune Infiltration Cell Type",decovolutionalMethod = "SVR", width =NA, height = NA, groupToCompare = NA, orderFactor = NA, colCount = 4){
 
   immuneMatrix =data.frame(immuneMatrix, Factor = groups)
    
  vars = colnames(immuneMatrix)[1:(dim(immuneMatrix)[2]-1)]
  a = list()
  testResults = data.frame()
  for  (i in 1:length(vars)){
    
    data.formula = formula(paste(vars[i] ,"~ Factor"))
    testResults.temp =  compare_means(data.formula,comparisons = project.groupsToCompare, p.adjust.method = "fdr", method='wilcox.test',data = immuneMatrix)
    testResults = rbind(testResults,testResults.temp)
    
    #  pairwise.wilcox.test(x = immuneMatrix[,i], g =  immuneMatrix$Factor, p.adjust.method = "fdr")
    
    #esto muestra pvalue y no padjusted es una cagada
    p2 = ggboxplot(immuneMatrix,  x ="Factor", y = vars[i],
                   title = gsub("\\."," ",vars[i]),
                   subtitle = " ",
                   fill = "Factor", add = c("none"),width=0.5,
                   order = orderFactor,
                   xlab = "",ylab ="TIL (%)",
                   repel = TRUE,combine = FALSE) +
      
      
      stat_compare_means(comparisons = groupToCompare,
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
  
  HELPER_SAVE_DATA_FRAME(as.data.frame(testResults),paste0(outputDir,"/testResults.csv" ))
  
  
  p3  = p3 +rremove("x.title") +rremove("x.axis")+rremove("y.text") +rremove("y.ticks")+rremove("y.axis")+rremove("y.title") 
  p3$labels$title = ""
  
  p3 = ggpar(p3,
             legend = "right", legend.title = "Group", font.legend =c(12,"bold", "black"))
  a[[i+1]] = p3
  
  b =   plot_grid(plotlist =a,
                  ncol = colCount)
  
  ggsave(b,file=paste0(outputDir,"SVR_BoxplotWithPValues.pdf"), width = width, height =height, limitsize = FALSE)

  
}
GET_IMMUNE_CIBERSORT_CLINICAL = function(immuneResults,samplesList){
  htseq.survival.fitlered = immuneResults
  rownames(htseq.survival.fitlered) = gsub("\\.","-",rownames(htseq.survival.fitlered))
  
  samples.clinical =  samplesList[samplesList$GROUP == "Primary Tumor",]#filtramos solo tumores
  
  htseq.survival.fitlered = htseq.survival.fitlered[rownames(htseq.survival.fitlered) %in% samplesList$Sample_ID,] 
  htseq.survival.fitlered = htseq.survival.fitlered[order(match (rownames(htseq.survival.fitlered), samplesList$Sample_ID)),] #reordenamos por las dudas
  
  samples.clinical = samplesList[samplesList$Sample_ID %in% rownames(htseq.survival.fitlered),]
  dim(samples.clinical)
  
  survivalData = GET_SURVIVAL_TIME(timeToLastFollow = samples.clinical$days_to_last_follow_up,
                                             vitalStatus = samples.clinical$vital_status, 
                                             daysToDeath = samples.clinical$days_to_death,
                                             otherVars = htseq.survival.fitlered)
  
  colnames(survivalData) = gsub(" ","_", colnames(survivalData))
  colnames(survivalData) = gsub("-","_", colnames(survivalData))
  colnames(survivalData) = gsub("\\+","", colnames(survivalData))
  colnames(survivalData) = gsub("\\(","", colnames(survivalData))
  colnames(survivalData) = gsub("\\)","", colnames(survivalData))
  
  return (survivalData)
}














