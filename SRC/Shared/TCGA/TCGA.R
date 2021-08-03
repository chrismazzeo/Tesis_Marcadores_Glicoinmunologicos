
#*******************************************************************************************************
#GENE LEVEL
#*******************************************************************************************************
# HTSeq + DESeq2
#*******************************************************************************************************

outputDir.htseq =  paste0(project.resultsPath, "HTSEQ/")
CreateGroupFolder(outputDir =outputDir.htseq,project.groupsToCompare)

 
htseq.rawCountMatrix = DESeq2_Merge_TCGA_ZippedExpressionFiles(sourceFolder = paste0(project.resourcePath,"rawData/htseq"),
                                                               sampleFolder = samplesList$File_ID,
                                                               samplesFileName =samplesList$File_Name,
                                                               sampleID = project.sample.name,
                                                               outputPath = NULL)

htseq.rawCountMatrix.annotation =  HELPER_RowNamesAsFirstColumn(htseq.rawCountMatrix, "ENSENMBL_GENE_ID")
htseq.rawCountMatrix.annotation = ANN_MergeAnnotationToDataFrame(htseq.rawCountMatrix.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")

HELPER_SAVE_DATA_FRAME(htseq.rawCountMatrix.annotation, paste0(project.resultsPath, "HTSEQ/CountMatrix/rawCountMatrixAnnotated.csv"))


htseq.dds = DESeq2_Matrix_Import(htseq.rawCountMatrix, samplesList, project.formula, coresNumber = setttings.coresNumber)
htseq.vsd = DESeq2_VST(htseq.dds)
HELPER_SAVE_DATA_FRAME(as.data.frame(assay(htseq.vsd)), paste0(project.resultsPath, "HTSEQ/CountMatrix/vsd.csv"))



if (project.rnaseq.run.da){
  DA_RNASEQ_ALL(rawCountMatrix = htseq.rawCountMatrix,
                deseq.dds = htseq.dds,
                dseq.vsd = htseq.vsd,
                samplesNames = project.sample.name,
                samplesGroups = project.sample.group,
                tsne.perplexity = 2,
                KMeans.kmax = length(levels(project.sample.group))*3,
                KMeans.optimalNumberOfClusters = length(levels(project.sample.group)),
                KMeans.maxIteration = settings.kmeans.numberOfinteractions,
                KMeans.numbersOfRandomCenters = settings.kmeans.numbersOfRandomCenters,
                groups.colorArray = settings.color,
                outputDir = outputDir.htseq)
}
#****************************************************************
#survival
if (project.rnaseq.run.expressionSurvival){
    htseq.survival.fitlered = assay(htseq.vsd)[ rownames(htseq.vsd) %in% geneFilteringList$ensembl_gene_id,]  
    htseq.survival.fitlered = t(htseq.survival.fitlered)
    
    samples.clinical =  samplesList[samplesList$GROUP == "Primary Tumor",]#filtramos solo tumores
    
    htseq.survival.fitlered = htseq.survival.fitlered[rownames(htseq.survival.fitlered) %in% samples.clinical$Sample_ID,] 
    htseq.survival.fitlered = htseq.survival.fitlered[order(match (rownames(htseq.survival.fitlered), samples.clinical$Sample_ID)),] #reordenamos por las dudas
    
    htseq.survivalData = GET_SURVIVAL_TIME(timeToLastFollow = samples.clinical$days_to_last_follow_up,
                                           vitalStatus = samples.clinical$vital_status,
                                           daysToDeath= samples.clinical$days_to_death,
                                          
                                           otherVars =  htseq.survival.fitlered)
    
   #pegamos el nombre del gen
    htseq.survivalData.id = data.frame(ensembl_gene_id = colnames(htseq.survivalData)[3:dim(htseq.survivalData)[2]])
    genes.id = data.frame(ensembl_gene_id = as.character(geneFilteringList$ensembl_gene_id), external_gene_name = as.character(geneFilteringList$external_gene_name), stringsAsFactors = FALSE)
    genes.names = merge (htseq.survivalData.id, genes.id, by.x ="ensembl_gene_id", by.y = "ensembl_gene_id")
    colnames(htseq.survivalData)[3:dim(htseq.survivalData)[2]] = paste0(genes.names$external_gene_name,"_", genes.names$ensembl_gene_id)
    
    # univariado
    GET_KAPLAN_MEIER_PLOT(survivalData = htseq.survivalData,
                      outputDir = paste0(outputDir.htseq,"survival/kaplan-meier/"))
    
    #multivariado
    #guardar
    try({
        GET_COXPH_PLOT (survivalData = htseq.survivalData,
                    outputDir = paste0(outputDir.htseq,"survival/coxph/"))
     

      })
    
  
    #por grupo 
    # sample.group.filtered = droplevels(project.sample.group[samplesList$GROUP == "Primary Tumor"])
    # group = levels(sample.group.filtered)
    # 
    # if (length(group) >1){
    #       for (i in 1:length(group)){
    #         htseq.survivalData.filtered.level = htseq.survivalData[sample.group.filtered == group[i],]
    # 
    #         GET_SURVIVAL_PLOT(survivalData = htseq.survivalData.filtered.level,
    #                       outputDir = paste0(outputDir.htseq,"survival/",group[i],"/"))
    #       
    #       }
    # }
      
}


#****************************************************************

#compare groups
for (i in 1:length(project.groupsToCompare)){
  contrast =project.groupsToCompare[[i]]
  outputDirGroup = paste0(outputDir.htseq,"DEA/",project.groupsToCompare[[i]][2],"_VS_",project.groupsToCompare[[i]][3],"/")
  
  htseq.res = DESeq2_ExtractResults (htseq.dds,filterFun = ihw, coresNumber = setttings.coresNumber, pAdjustMethod = "BH", adjustedPValue = project.settings.DEA.padj,contrast,  format ="DataFrame", lfcThreshold = 0)
 # htseq.shrunken <- lfcShrink(htseq.dds,coef = resultsNames(htseq.dds)[2], lfcThreshold = project.settings.DEA.H0.Log2FC)  #SOLO para plotear y raankear, despues da lo mismo
  
  htseq.res.annotation = HELPER_RowNamesAsFirstColumn(htseq.res, "ENSENMBL_GENE_ID")
  htseq.res.annotation = ANN_MergeAnnotationToDataFrame(htseq.res.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")
  
 # hstseq.shrunken.annotation = HELPER_RowNamesAsFirstColumn(htseq.shrunken, "ENSENMBL_GENE_ID")
#  hstseq.shrunken.annotation = ANN_MergeAnnotationToDataFrame(hstseq.shrunken.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")
  
  # DA_PLOT_TCGA_GENE_EXPRESSION (expressionMatrix = htseq.vsd,
  #                            samplesList = samplesList,
  #                            geneFilteringList = geneFilteringList,
  #                            group = project.sample.group,
  #                            outputDir = paste0(outputDirGroup, "GeneExpression/")
  # )
  
  DESeq2_Result_ALL_Stats(DESeqResult = htseq.res, 
                          DESeqResult.shrunken =  NULL, 
                          DESeqResult.annotation = htseq.res.annotation,
                          padj = project.settings.DEA.padj,
                          log2FCCutOff = project.settings.DEA.cutoff.log2FC, 
                          outputDir = outputDirGroup
  )
  
  HELPER_FILTER_SAVE_PFI (DESeqResult.annotation = htseq.res.annotation,
                          geneList = geneFilteringList,
                          padj = project.settings.DEA.padj, 
                          log2FC = project.settings.DEA.cutoff.log2FC, 
                          outputDir = outputDirGroup 
  )
  #log2FC = 1
  if  (project.settings.DEA.H0.Log2FC >0){
        
        outputDirGroup = paste0(outputDir.htseq,"DEA/",project.groupsToCompare[[i]][2],"_VS_",project.groupsToCompare[[i]][3],"_Log2FC_",project.settings.DEA.H0.Log2FC,"/")
       
        htseq.res = DESeq2_ExtractResults (htseq.dds,filterFun = ihw, coresNumber = setttings.coresNumber, pAdjustMethod = "BH", adjustedPValue = project.settings.DEA.padj,contrast,  format ="DataFrame", lfcThreshold = project.settings.DEA.H0.Log2FC)
        # htseq.shrunken <- lfcShrink(htseq.dds,coef = resultsNames(htseq.dds)[2], lfcThreshold = project.settings.DEA.H0.Log2FC)  #SOLO para plotear y raankear, despues da lo mismo
        
        htseq.res.annotation = HELPER_RowNamesAsFirstColumn(htseq.res, "ENSENMBL_GENE_ID")
        htseq.res.annotation = ANN_MergeAnnotationToDataFrame(htseq.res.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")
        
        # hstseq.shrunken.annotation = HELPER_RowNamesAsFirstColumn(htseq.shrunken, "ENSENMBL_GENE_ID")
        #  hstseq.shrunken.annotation = ANN_MergeAnnotationToDataFrame(hstseq.shrunken.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")
        
    
        
        DESeq2_Result_ALL_Stats(DESeqResult = htseq.res, 
                                DESeqResult.shrunken =  NULL, 
                                DESeqResult.annotation = htseq.res.annotation,
                                padj = project.settings.DEA.padj,
                                log2FCCutOff = project.settings.DEA.cutoff.log2FC, 
                                outputDir = outputDirGroup
        )
      
       HELPER_FILTER_SAVE_PFI (DESeqResult.annotation = htseq.res.annotation,
                              geneList = geneFilteringList,
                              padj = project.settings.DEA.padj, 
                              log2FC = project.settings.DEA.cutoff.log2FC, 
                              outputDir = outputDirGroup 
      )
       
      
  }
}
















