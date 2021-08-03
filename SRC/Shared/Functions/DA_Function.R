#***********************************************
# SampleDistance
#***********************************************

DA_PlotSamplesDistance = function (expressionMatrix, sampleNames,sampleGroups, size = 10,filePath){
  
  sampleDistance <-  dist(expressionMatrix)
  sampleDistance = as.matrix(sampleDistance)
  rownames(sampleDistance) = sampleNames
#  colnames(sampleDistance) = NULL

 rowAnnotation <- data.frame(GROUP = sampleGroups)
  rownames(rowAnnotation) = sampleNames
 
  annotationCol =data.frame(Type = sampleGroups)
  rownames(annotationCol) = sampleNames
  
  type.colors = settings.graphics.colors[1:length(levels(factor(sampleGroups)))]
  names(type.colors) = levels(factor(sampleGroups))
  
  family.colors = settings.graphics.colors[3:(length(levels(sampleGroups))+2)]
  names(family.colors) = levels(sampleGroups)
  annotationColor = list(Type = type.colors,
                         Family = family.colors)
  
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  graphics =  pheatmap(sampleDistance,
                       cellwidth = size,cellheight = size,
                       main = "Sample Distance\n",
                       col=colors, 
                       annotation_col = annotationCol,annotation_colors = annotationColor,
                       show_rownames = FALSE,
                       treeheight_row = 0)
  if(!is.null(filePath)){
    pdf(filePath, width = 10)
    print(graphics)
    dev.off()
    
    print(graphics)
  }
}
#***********************************************
#NOTA = 1. hay que escalar los datos con vst rlog, log(n+1), scale 
#       2. las muestras tienen que estar como filas
#***********************************************
# TSNE
#***********************************************
# ploteo con tsne, agrego el label, y com kmeans coloreo los grupos  o coloreo segun el sample type y que quede como quede
DA_GetTSNE = function (dataFrame,dimensions, perplexity = 500, max_iter = 1000, checkDuplicates = TRUE){
  set.seed(123)
  tsne <- Rtsne(dataFrame, dims = dimensions, perplexity=perplexity, verbose=TRUE, max_iter = max_iter, check_duplicates = checkDuplicates)
  return (tsne)
}
#***********************************************
# KMEANS
#***********************************************
#type = wss ,silhouette,gap_stat,NbClust
DA_GetOptimalNumberClusters = function (dataFrame,type,kMax = 2, filePath = NA) {
  
  if (type == "wss" | type == "all"){ 
    optimalCluster = fviz_nbclust(dataFrame, kmeans, method = "wss",k.max = kMax) +
      labs(subtitle = "Elbow method")
    print(optimalCluster)
  }
  
  if (type == "silhouette" | type == "all"){
    optimalCluster = fviz_nbclust(dataFrame, kmeans, method = "silhouette",k.max = kMax) +
      labs(subtitle = "Silhouette method")
    print(optimalCluster)
  }
  
  if (type == "gap_stat" | type == "all"){
    optimalCluster = fviz_nbclust(dataFrame, kmeans, method = "gap_stat",nboot = 50,k.max = kMax) +
      labs(subtitle = "Gap statistic method")
    print(optimalCluster)
  }
  
  if (type == "NbClust" | type == "all"){
    nbClust <- NbClust(dataFrame, distance = "euclidean", min.nc = 1,max.nc = kMax, method = "kmeans")
    optimalCluster = fviz_nbclust(nbClust)+ labs(subtitle = "NbClust")
    print(optimalCluster)
  }
  
  if (!is.na(filePath)){
    HELPER_SAVE_PDF(optimalCluster, filePath)
  }
  return (optimalCluster)
}

DA_GetKmeans = function (dataFrame, numberOfCenters = 2, maxIteration= 100, numberOfRandomCenters = 1){
  clusters<-kmeans(dataFrame, centers=numberOfCenters,iter.max = maxIteration, nstart = numberOfRandomCenters)  
  print( table(clusters$cluster))
  return (clusters)
}


#***********************************************
# PCA
#***********************************************
DA_GetPCA = function(dataFrame, scale =FALSE){
  
  pca = PCA(dataFrame, scale.unit = scale, ncp = 5, graph = FALSE)
  summary(pca)
  fviz_eig(pca, addlabels=TRUE, hjust = -0.3, title =  "Percentage explained by each PC")+ theme_minimal()
  return (pca)
}
DA_PlotPCA = function(pca){
  fviz_pca_ind(pca)
}

DA_PlotSamplesCoockDistance = function (dds,samplesNames,outputDir =  NULL){
  
  #outliers,
  #cuando hay >=3 muestras, se calcula la distancia de cook`s, si es mayor a cierto treshold se pone como outlier, y el pvalue y padj es NA`
  #With many degrees of freedom – i.,e., many more samples than number of parameters to be estimated – it is undesirable to remove entire genes from the analysis just because their data include a single count outlier. When there are 7 or more replicates for a given sample, the DESeq function will automatically replace counts with large Cook’s distance with the trimmed mean over all samples, scaled up by the size factor or normalization factor for that sample. This approach is conservative, it will not lead to false positives, as it replaces the outlier value with the value predicted by the null hypothesis. This outlier replacement only occurs when there are 7 or more replicates, and can be turned off with DESeq(dds, minReplicatesForReplace=Inf).
  #summary(res) da cuantos oluliers
  
  #si da muchos, hay que analizar si hay alguna muestra  es un outlier
  #When there are thousands of reported outliers, it might make more sense to turn off the outlier filtering/replacement (DESeq with minReplicatesForReplace=Inf and results with cooksCutoff=FALSE) and perform manual inspection: First it would be advantageous to make a PCA plot as described above to spot individual sample outliers; Second, one can make a boxplot of the Cook’s distances to see if one sample is consistently higher than others (here this is not the case):
  pdf(paste0(outputDir,"CooksDistance.pdf"))
  par(mar=c(8,5,2,2))
  dataframe = log10(assays(dds)[["cooks"]])
  colnames(dataframe) = samplesNames
  
  boxplot(dataframe, range=0, las=2, main = "Cook`s distance")
  dev.off()
  
  pdf(paste0(outputDir,"Disperssion.pdf"))
  plotDispEsts(dds)
  dev.off()
  
  
}
DA_HEATMAP_CountMatrix = function(dds,vsd,sampleNames, sampleGroups,filePath){
  
  #sort by expression
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)
  
  df <- data.frame(GROUP = sampleGroups)
  rownames(df) = sampleNames
  colnames(vsd) = sampleNames
  graphics = pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                      cluster_cols=TRUE, annotation_col = df)
  
  
  HELPER_SAVE_PDF(graphics,filePath)
  print(graphics)
  
}
DA_MICROARRAY_HEATMAP_CountMatrix = function(rawExpressionMatrix,sampleNames, sampleGroups,filePath){
  
  #sort by expression
  select <- order(rowMeans(log(rawExpressionMatrix,2)),
                  decreasing=TRUE)
  
  df <- data.frame(GROUP = sampleGroups)
  rownames(df) = sampleNames
  
  colnames(rawExpressionMatrix) = sampleNames
  rawExpressionMatrix = rawExpressionMatrix[select,]
  graphics = pheatmap(rawExpressionMatrix, cluster_rows=FALSE, show_rownames=FALSE,
                      cluster_cols=TRUE, annotation_col = df)
  
  
  HELPER_SAVE_PDF(graphics,filePath)
  print(graphics)
  
}

#heatmap of count matrix
#***********************************************
# Plots
#***********************************************
DA_Plot2D = function(dataFrame, title ="", xTitle ="",yLabel ="", colorsArray = NA, group = NA, outputDir, fileName){
  
  if (!is.na(group)){
    factorColor = HELPER_ConvertFactorToInteger(group)
    colorsArray = colorsArray[factorColor]
  }
  else{
    if(is.na(colorsArray)){
      colorsArray =  "black"
    }
  }
  
  
  
  
  graphics <- plot_ly(x = dataFrame[,1], y = dataFrame[,2], color = colorsArray, name =  factor(group), text = group)%>%
    layout(title = title,
           xaxis = list(title = xTitle,
                        zeroline = TRUE),
           
           yaxis = list(title =yLabel)
    )
  HELPER_SAVE_ORCA(graphics, outputDir = outputDir, fileName = fileName)
  
  return (graphics)
  
}
DA_Plot3D = function(dataFrame, title ="", xTitle ="",yLabel ="", colorsArray = NA, group = NA, outputDir, fileName){
  
  if (!is.na(group)){
    factorColor = HELPER_ConvertFactorToInteger(group)
    colorsArray = colorsArray[factorColor]
  }
  else{
    if(is.na(colorsArray)){
      colorsArray =  "black"
    }
  }
  
  graphics <- plot_ly(x = dataFrame[,1], y = dataFrame[,2],z = dataFrame[,3], color = colorsArray, name =  factor(group), text = group)%>%
    layout(title = title,
           xaxis = list(title = xTitle,
                        zeroline = TRUE),
           
           yaxis = list(title =yLabel)
    )
  
  
  HELPER_SAVE_ORCA(graphics, outputDir = outputDir, fileName = fileName)
  return (graphics)
}

#*******************************************************************************************************
#KMeans.kmax depende de la cantidad de muestras que halla

DA_RNASEQ_ALL = function (rawCountMatrix, deseq.dds,dseq.vsd, samplesNames, samplesGroups, tsne.perplexity,KMeans.kmax,KMeans.optimalNumberOfClusters,KMeans.maxIteration = 100,KMeans.numbersOfRandomCenters, groups.colorArray,outputDir){
  
  expressionMatrix = t(assay(dseq.vsd))
  rownames(expressionMatrix) = samplesNames
  
  #******************************************************************
  #sample-to-sample distance  muestra que samples estan mas cercanas
  #*****************************************************************
  rawCountMatrix2 = t(rawCountMatrix)
  rownames(rawCountMatrix2) = samplesNames
  DA_PlotSamplesDistance(rawCountMatrix2,sampleNames = samplesNames, sampleGroups= samplesGroups, filePath = paste0(outputDir,"DescriptiveAnalysis/sampleDistance.pdf"))
  
  #******************************************************************
  #Cook´s distance
  #******************************************************************
  DA_PlotSamplesCoockDistance(deseq.dds,samplesNames,outputDir = paste0(outputDir,"DescriptiveAnalysis/"))
  #******************************************************************
  #vemos raw expresion
  #*****************************************************************
  DA_HEATMAP_CountMatrix(deseq.dds,dseq.vsd, samplesNames, samplesGroups,filePath =  paste0(outputDir,"DescriptiveAnalysis/vsdHeatmap.pdf"))
  
  #**************
  #KMEANS
  #**************
  optimalCluster.elbow = DA_GetOptimalNumberClusters(dataFrame = expressionMatrix,type = "wss",kMax = KMeans.kmax, filePath = paste0(outputDir,"DescriptiveAnalysis/Elbow.pdf"))
  optimalCluster.silhoutte = DA_GetOptimalNumberClusters(dataFrame = expressionMatrix,type = "silhouette",kMax = KMeans.kmax, filePath = paste0(outputDir,"DescriptiveAnalysis/Silhouette.pdf"))
  optimalCluster.gap_stat = DA_GetOptimalNumberClusters(dataFrame = expressionMatrix,type = "gap_stat",kMax = KMeans.kmax, filePath = paste0(outputDir,"DescriptiveAnalysis/GapSatistics.pdf"))
  optimalCluster.nbClust = DA_GetOptimalNumberClusters(dataFrame = expressionMatrix,type = "NbClust",kMax = 3) #se cuelga
  
   clusters.kmeans = DA_GetKmeans(expressionMatrix, numberOfCenters = KMeans.optimalNumberOfClusters, maxIteration = KMeans.maxIteration, numberOfRandomCenters = KMeans.numbersOfRandomCenters)
  #save clusters
   samplesList = data.frame(samplesNames = samplesNames, KMEANS = clusters.kmeans$cluster)
   HELPER_SAVE_DATA_FRAME(samplesList, paste0(outputDir,"DescriptiveAnalysis/KMeansClusters.csv")) #raw htseq matrix
  #**************
  #TSNE
  #**************
  tsne.2d = DA_GetTSNE(dataFrame = expressionMatrix, dimensions = 2,perplexity = tsne.perplexity)
  #clusters.tsne2d = DA_GetKmeans(tsne.2d$Y, numberOfCenters = KMeans.optimalNumberOfClusters, maxIteration = KMeans.maxIteration, numberOfRandomCenters = KMeans.numbersOfRandomCenters)
  DA_Plot2D(tsne.2d$Y, title = paste0("TSNE 2D - Clustered "), colorsArray = groups.colorArray, group = samplesGroups,outputDir = paste0(outputDir,"DescriptiveAnalysis/"), fileName = "TSNE2d.pdf")
  
  tsne.3d = DA_GetTSNE(dataFrame = expressionMatrix, dimensions = 3,perplexity = tsne.perplexity)
  #clusters.tsne3d = DA_GetKmeans(tsne.3d$Y, numberOfCenters = KMeans.optimalNumberOfClusters, maxIteration = KMeans.maxIteration, numberOfRandomCenters = KMeans.numbersOfRandomCenters)
  DA_Plot3D(tsne.3d$Y, title = "TSNE 3D -Clustered", colorsArray = groups.colorArray, group = samplesGroups, outputDir = paste0(outputDir,"DescriptiveAnalysis/"), fileName = "TSNE3d.pdf")
  scatter3d(tsne.3d$Y[,1],tsne.3d$Y[,2], tsne.3d$Y[,3],point.col = "blue", surface=FALSE, groups = samplesGroups,ellipsoid = TRUE, grid = FALSE,axis.scales = FALSE,xlab = "", ylab = "", zlab = "")
  
  rgl.postscript( paste0(outputDir,"DescriptiveAnalysis/tsne3d_2"),fmt="pdf")
  #**************
  #PCA
  #**************
  
  plotPCA(dseq.vsd, intgroup=c("GROUP")) #deseq2PCA
  
  pca = DA_GetPCA(expressionMatrix)
  DA_PlotPCA(pca)
}

DA_MICROARRAY_ALL = function (rawCountMatrix,samplesNames, samplesGroups, tsne.perplexity,KMeans.kmax,KMeans.optimalNumberOfClusters,KMeans.maxIteration = 100,KMeans.numbersOfRandomCenters, groups.colorArray,outputDir){
  
 
  expressionMatrix = t(exprs(rawCountMatrix))
  expressionMatrix = log(expressionMatrix,2)
  rownames(expressionMatrix) = samplesNames
  
  #******************************************************************
  #sample-to-sample distance  muestra que samples estan mas cercanas
  #*****************************************************************
  DA_PlotSamplesDistance(expressionMatrix,sampleNames = samplesNames, sampleGroups= samplesGroups, filePath = paste0(outputDir,"DescriptiveAnalysis/sampleDistance.pdf"))
  
  #******************************************************************
  #vemos raw expresion
  #*****************************************************************
  DA_MICROARRAY_HEATMAP_CountMatrix(exprs(rawCountMatrix), samplesNames, samplesGroups,filePath =  paste0(outputDir,"DescriptiveAnalysis/ExpressionHeatmap.pdf"))
  
  #**************
  #KMEANS
  #**************
  optimalCluster.elbow = DA_GetOptimalNumberClusters(dataFrame = expressionMatrix,type = "wss",kMax = KMeans.kmax, filePath = paste0(outputDir,"DescriptiveAnalysis/Elbow.pdf"))
  optimalCluster.silhoutte = DA_GetOptimalNumberClusters(dataFrame = expressionMatrix,type = "silhouette",kMax = KMeans.kmax, filePath = paste0(outputDir,"DescriptiveAnalysis/Silhouette.pdf"))
  optimalCluster.gap_stat = DA_GetOptimalNumberClusters(dataFrame = expressionMatrix,type = "gap_stat",kMax = KMeans.kmax, filePath = paste0(outputDir,"DescriptiveAnalysis/GapSatistics.pdf"))
  #optimalCluster.nbClust = DA_GetOptimalNumberClusters(dataFrame = expressionMatrix,type = "NbClust",kMax = 3) #se cuelga
  
  clusters.kmeans = DA_GetKmeans(expressionMatrix, numberOfCenters = KMeans.optimalNumberOfClusters, maxIteration = KMeans.maxIteration, numberOfRandomCenters = KMeans.numbersOfRandomCenters)
  #save clusters
  samplesList = data.frame(samplesNames = samplesNames, KMEANS = clusters.kmeans$cluster)
  HELPER_SAVE_DATA_FRAME(samplesList, paste0(outputDir,"DescriptiveAnalysis/KMeansClusters.csv"), rowNames = FALSE) #raw htseq matrix
  #**************
  #TSNE
  #**************
  tsne.2d = DA_GetTSNE(dataFrame = expressionMatrix, dimensions = 2,perplexity = tsne.perplexity)
  #clusters.tsne2d = DA_GetKmeans(tsne.2d$Y, numberOfCenters = KMeans.optimalNumberOfClusters, maxIteration = KMeans.maxIteration, numberOfRandomCenters = KMeans.numbersOfRandomCenters)
  DA_Plot2D(tsne.2d$Y, title = paste0("TSNE 2D - Clustered "), colorsArray = groups.colorArray, group = samplesGroups,outputDir = paste0(outputDir,"DescriptiveAnalysis/"), fileName = "TSNE2d.pdf")
  
  tsne.3d = DA_GetTSNE(dataFrame = expressionMatrix, dimensions = 3,perplexity = tsne.perplexity)
  #clusters.tsne3d = DA_GetKmeans(tsne.3d$Y, numberOfCenters = KMeans.optimalNumberOfClusters, maxIteration = KMeans.maxIteration, numberOfRandomCenters = KMeans.numbersOfRandomCenters)
  DA_Plot3D(tsne.3d$Y, title = "TSNE 3D -Clustered", colorsArray = groups.colorArray, group = samplesGroups, outputDir = paste0(outputDir,"DescriptiveAnalysis/"), fileName = "TSNE3d.pdf")
  
  #**************
  #PCA
  #**************
  
  pca = DA_GetPCA(expressionMatrix)
  DA_PlotPCA(pca)
}


DA_PLOT_TCGA_GENE_EXPRESSION = function(expressionMatrix, geneFilteringList, samplesList, group,  outputDir){
  
  dir.create(file.path( outputDir), recursive = TRUE, showWarnings = FALSE)
  expressionMatrix.filtered = as.data.frame(assay(expressionMatrix))
  expressionMatrix.filtered = expressionMatrix.filtered[rownames(expressionMatrix) %in% geneFilteringList$ensembl_gene_id,]
  dim(expressionMatrix.filtered)
  #htseq.vsd.filtered = t(assay(htseq.vsd.filtered))
  expressionMatrix.filtered = as.data.frame(t(expressionMatrix.filtered))
  
  #***** heatmap
  rowAnnotation2 <- data.frame(Type = samplesList$GROUP, Stage = samplesList$tumor_stage, MSI = samplesList$msi_status, Lost_Repair = samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc)
  rownames(rowAnnotation2) = rownames(expressionMatrix.filtered)
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
  
  
  pheatmap(t(expressionMatrix.filtered), main = "Glyco Genes expression", cluster_rows = FALSE, cluster_cols = FALSE, 
           annotation_col = rowAnnotation2 ,  treeheight_col = 0,show_rownames = TRUE,show_colnames = FALSE,
           annotation_colors = ann_colors,
           cellwidth = 3,cellheight = 3, fontsize_col = 3,fontsize_row = 3,
           
           filename =paste0(outputDir, "geneExpressionHeatmap.pdf"))
  #************************ boxplot 
  #expressionMatrix.filtered = as.data.frame(t(expressionMatrix.filtered))
  
  # expressionMatrix.filtered$Group = samplesList$GROUP
  expressionMatrix.filtered$Group = group
  expressionMatrix.filtered.melt = melt(expressionMatrix.filtered,   id.vars = "Group")
  
  family = levels(geneFilteringList$Group)
  outputDir = paste0(outputDir,"boxplot/")
  dir.create(file.path( outputDir), recursive = TRUE, showWarnings = FALSE)
  for (i in 1:length(family)){
    expressionMatrix.filtered.melt.temp = expressionMatrix.filtered.melt[expressionMatrix.filtered.melt$variable %in% geneFilteringList$ensembl_gene_id[geneFilteringList$Group == family[i]],]
    
    p <- ggplot(data = expressionMatrix.filtered.melt, aes(x=variable, y=value)) 
    p <- p + geom_boxplot(aes(fill=Group))
    
    p <- p + facet_wrap( ~ variable, scales="free")
    p <- p + xlab("") + ylab("") + ggtitle(family[i])
    p <- p + guides(fill=guide_legend(title="Type"))
    p =  p +settings.graphics.theme +settings.graphics.borderColor
    fileName = gsub(" ","_",family[i])
    fileName = gsub("/","_",fileName)
    ggsave(paste0(outputDir,"./",fileName,".pdf"),p )
  }
}
