

#****************************
#load/save data frames to TSV
#****************************
HELPER_LoadTSV = function(inputFile){
  samples = read_tsv(inputFile)
  samples = as.data.frame(samples)
  return (samples)
}
HELPER_SaveTSV = function (dataFrame, outputPath){
  if (!is.null(outputPath)){
    write_tsv(dataFrame,path = outputPath )
  }
}

HELPER_LoadCSV = function(inputFile, colNames = TRUE){
  samples = read.csv(inputFile, header = colNames)
  samples = as.data.frame(samples)
  return (samples)
}
HELPER_SaveCSV = function(dataFrame, outputPath, rowNames = TRUE){
  if (!is.null(outputPath)){

      write.csv(dataFrame, outputPath, row.names=rowNames, quote=F)
  }
}

HELPER_SAVE_EXCEL= function(dataFrame, outputPath, rowNames = TRUE){ 

  if (!is.null(outputPath)){
   
   
    write_xlsx(dataFrame, outputPath, col_names = TRUE, format_headers = TRUE)
  }
}
HELPER_SAVE_DATA_FRAME = function (dataFrame, outputPath, rowNames = TRUE, type = "TSV"){
  if (!is.null(outputPath)){
      if (type == "TSV"){
        HELPER_SaveTSV(dataFrame, outputPath)
      }
      if (type == "CSV"){
        HELPER_SaveCSV(dataFrame, outputPath, rowNames)
      }
      
      outputPath  = gsub("csv", "xlsx",outputPath)
      outputPath  = gsub("tsv", "xlsx",outputPath)
      HELPER_SAVE_EXCEL(dataFrame, outputPath, rowNames)
  }
}


HELPER_InsertAtFirstColumn = function(dataFrame, column, firstColumName){
  dataFrame = data.frame(fistColumnName = column, dataFrame)
  return (dataFrame)
}
HELPER_FirstColumnAsRowName = function(dataFrame){
  rownames(dataFrame) = dataFrame[,1]
  dataFrame = dataFrame[,-1]
  return (dataFrame)
  
}
HELPER_RowNamesAsFirstColumn = function(dataFrame, fistColumnName){
  dataFrame = data.frame(a = rownames(dataFrame), dataFrame)
  colnames(dataFrame)[1] = fistColumnName
  
  return (dataFrame)
}
HELPER_ConvertFactorToInteger = function (factor){
  sampleFactor = as.character(factor)

  sampleLevels = levels(as.factor(factor))
  for (i in 1:length(sampleLevels)){
    sampleFactor[sampleFactor == sampleLevels[i]] = i
  }
  sampleFactor = as.numeric(sampleFactor)
  return (sampleFactor)
}

HELPER_Sort = function(dataFrame, rowOrder, decreasing=FALSE){
  
  return( dataFrame[order(rowOrder,decreasing = decreasing),])
  
}

HELPER_SAVE_PDF = function(graphic, filePath, draw = TRUE){

  pdf(filePath)
  print(graphic)
  dev.off()
}

HELPER_FILTER_SAVE_PFI = function( DESeqResult.annotation,geneList, padj = 0.05, log2FC = 0, outputDir, project.type = "RNASEQ"){
  dir.create(file.path(paste0(outputDir,"DEA/GlycoHeatmaps/Log2FC_1")), recursive = TRUE, showWarnings = FALSE)
  
  HELPER_SAVE_DATA_FRAME(DESeqResult.annotation, paste0(outputDir,"DEA/RES_ALL_GENES.csv"), rowNames = FALSE) #Todos
  
  res.filtered.significant = FILTERING_DESEQ2_RES_ByAdjPvalue(DESeqResult.annotation,padj)#solo significativos
  HELPER_SAVE_DATA_FRAME(res.filtered.significant,  paste0(outputDir,"DEA/RES_SIGNIFICANT_GENES.csv"), rowNames = FALSE)
  
  if(log2FC != 0){  # significativos y |log2FC| >1
    res.filtered.significant.log2FC = res.filtered.significant[abs(res.filtered.significant$log2FoldChange) > log2FC,]  
    HELPER_SAVE_DATA_FRAME(res.filtered.significant.log2FC,paste0(outputDir,"DEA/RES_SIGNIFICANT_GENES_LOG2FC.csv"), rowNames = FALSE) 
  }
 
 ######## GeneList ######
  if (project.type == "RNASEQ"){
    res.filtered.allglyco = FILTERING_ByGeneList(DESeqResult.annotation,DESeqResult.annotation$ENSENMBL_GENE_ID,geneList$ensembl_gene_id) #todos glyco
    #agregamos la familia
    tempData = data.frame (ensembl_gene_id = geneList$ensembl_gene_id, Family = geneList$Group)
    res.filtered.allglyco = merge(tempData, res.filtered.allglyco, by.x = "ensembl_gene_id" , by.y = "ENSENMBL_GENE_ID")
  }
  else{
    res.filtered.allglyco = FILTERING_ByGeneList(DESeqResult.annotation,DESeqResult.annotation$external_gene_name,geneList$external_gene_name) #todos glyco
    
    tempData = data.frame (external_gene_name = geneList$external_gene_name, Family = geneList$Group)
    res.filtered.allglyco = merge(tempData, res.filtered.allglyco, by.x ="external_gene_name" , by.y = "external_gene_name", all.x = TRUE )
  }

  HELPER_SAVE_DATA_FRAME(res.filtered.allglyco,paste0(outputDir,"DEA/RES_ALL_GLYCO_GENES.csv"), rowNames = FALSE)

  
  res.filtered.significant.glyco = FILTERING_DESEQ2_RES_ByAdjPvalue(res.filtered.allglyco,padj) #glyco significativos
  HELPER_SAVE_DATA_FRAME(res.filtered.significant.glyco,paste0(outputDir,"DEA/RES_SIGNIFICANT_GLYCO_GENES.csv"), rowNames = FALSE)
  
  if(log2FC != 0){ #glyco significativos y |log2FC| >1 
    res.filtered.significant.glyco.log2FC = res.filtered.significant.glyco[abs(res.filtered.significant.glyco$log2FoldChange) > log2FC,] 
    HELPER_SAVE_DATA_FRAME(res.filtered.significant.glyco.log2FC,paste0(outputDir,"DEA/RES_SIGNIFICANT_GLYCO_GENES_LOG2FC.csv"), rowNames = FALSE)
  }
  #*******************************************************************************************************
  #Glyco heatmap
  #*******************************************************************************************************
  if (dim(res.filtered.significant.glyco)[1] > 0){
    HELPER_GENES_HEATMAP (res.filtered.significant.glyco, paste0(outputDir,"DEA/GlycoHeatmaps/"))
  }
  if (dim(res.filtered.significant.glyco.log2FC)[1] > 0){
    HELPER_GENES_HEATMAP (res.filtered.significant.glyco.log2FC, paste0(outputDir,"DEA/GlycoHeatmaps/Log2FC_1/"))
  }
  
  
  return (list(significant.all = res.filtered.significant,
               significant.all.log2FC = res.filtered.significant.log2FC,
               filtered = res.filtered.allglyco ,
               filtered.significant = res.filtered.significant.glyco,
               filtered.significant.log2FC = res.filtered.significant.glyco.log2FC))
}
HELPER_GENES_HEATMAP = function (glycoGenes, outputDir){
  
  tempGlyco =glycoGenes
  tempGlyco = tempGlyco[order(tempGlyco$log2FoldChange, decreasing = TRUE),]
  
  tempGlyco.all = as.matrix(tempGlyco$log2FoldChange)
  rownames(tempGlyco.all) = tempGlyco$external_gene_name
  
  #all
  try(
     pheatmap(tempGlyco.all, 
           cluster_rows = FALSE,cluster_cols = FALSE, cellwidth = 30, cellheight =10,display_numbers = TRUE,na_col = "black",
           filename = paste0(outputDir,"ALL_GLYCO_GENES_HEATMAP.pdf"))
  
  )
  breaksList = seq(-8, 8, by = 2)
  
  #per family
  for (i in 1:length(levels(factor(tempGlyco$Family)))){
    tempGlyco.family = tempGlyco[tempGlyco$Family == levels(factor(tempGlyco$Family))[i],]
    tempGlyco.all = as.matrix(tempGlyco.family$log2FoldChange)
    rownames(tempGlyco.all) = tempGlyco.family$external_gene_name
    if (dim(tempGlyco.all)[1]>0){ #no dibujamos las familias que no tienen ninguno
    
        
        title = levels(factor(tempGlyco$Family))[i]
        
        
        title = gsub("/", " ", title)
      
        cellheight =10
        try(
          
          pheatmap(tempGlyco.all, 
                   
                   fontsize =5,
                   color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaksList)),
                   breaks = breaksList,
                   cluster_rows = FALSE,cluster_cols = FALSE, cellwidth = 30, cellheight =cellheight,display_numbers = TRUE, 
                    width = 6, height = length(tempGlyco.all)/6+1,
                   na_col = "black",
                   filename = paste0(outputDir, title,".pdf"))
        )
    
    }
    
  }
  
  try({
      tempGlyco =glycoGenes
      tempGlyco = tempGlyco[order(tempGlyco$log2FoldChange, decreasing = TRUE),]
      
      tempGlyco.all = as.matrix(tempGlyco$log2FoldChange)
      rownames(tempGlyco.all) = tempGlyco$external_gene_name
      grid.arrange(rectGrob(), rectGrob())
      max.row = 50
      max.col = 3
      max.pages = ceiling(length(tempGlyco.all)/max.row/max.col)
      
      col.number = ceiling(length(tempGlyco.all)/max.row)
      
      
      for (pageNumber in 1:max.pages){
        plot_list=list()
        start =  max.row*max.col*(pageNumber-1)+pageNumber
        offset = start +max.row*max.col
        if (offset> length(tempGlyco.all)){
          offset = length(tempGlyco.all) 
        }
        
        col.number = ceiling((offset-start)/max.row)
        for (i in 1:col.number){
          pos.start = start + max.row*(i-1)
          if ((pos.start + max.row) > length(tempGlyco.all)){
            post.end = length(tempGlyco.all)
          }else{
            post.end = pos.start + max.row
          }
          
          temp = as.matrix( tempGlyco.all[pos.start:post.end,])
          
          
          x = pheatmap(temp, 
                       color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaksList)),
                       
                       breaks = breaksList,
                       cluster_rows = FALSE,cluster_cols = FALSE, cellwidth = 40, cellheight =14,display_numbers = TRUE,fontsize_number = 10,fontsize_row = 13,
                       na_col = "black", silent = TRUE)
          
          plot_list[[i]] = x[[4]]
          
        }
        
        g <- grid.arrange(arrangeGrob(grobs= plot_list,ncol=col.number))
        ggsave(paste0(outputDir,"ALL_GLYCO_GENES_HEATMAP_",pageNumber,".pdf"), g, height = 11)
 
    
    }    
  })

  
  # plot_list=list()
  #for (i in 1:length(levels(tempGlyco$Family))){
  # tempGlyco.family = tempGlyco[tempGlyco$Family == levels(tempGlyco$Family)[i],]
  #tempGlyco.all = as.matrix(tempGlyco.family$log2FoldChange)
  # rownames(tempGlyco.all) = tempGlyco.family$external_gene_name
  
  
  # title = levels(tempGlyco$Family)[i]
  
  
  # title = gsub("/", " ", title)
  # print (title)
  
  #  x= pheatmap(tempGlyco.all,cluster_rows = FALSE,cluster_cols = FALSE, cellwidth = 30, cellheight =10,display_numbers = TRUE,
  #             main = title,
  #            fontsize =5,
  #           colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(breaksList)),
  #          breaks = breaksList
  #         )
  #plot_list[[i]] = x[[4]]
  
  
  #}
  #g <- grid.arrange(arrangeGrob(grobs= plot_list,ncol=19))
  #ggsave("./g.pdf",g, width = 40, height = 30)
  #HELPER_SAVE_PDF(glycoGraphics, paste0(settings.folder.output.dea,project.program,"_SIGNIFICANT_GLYCO_GENES_HEATMAP.PDF") )
  
}

HELPER_FILER_SAVE_DEXEQ_PFI = function (resDataFrame,geneList, padj = 0.05, outputDir){
  
  resDataFrame = as.data.frame(resDataFrame)
  resDataFrame = HELPER_RowNamesAsFirstColumn(resDataFrame,"ID")
  resDataFrame$transcripts = as.character(resDataFrame$transcripts)
  resDataFrame$transcripts = gsub(","," - ",resDataFrame$transcripts)
  
  HELPER_SAVE_DATA_FRAME(as.data.frame(resDataFrame),paste0(outputDir,"/DEA/RES_ALL_GENES.csv") )
  
  
  dexseq.res.significantID = unique(resDataFrame$groupID[!is.na(resDataFrame$padj) &resDataFrame$padj <padj])
  dexseq.res.significant = resDataFrame[resDataFrame$groupID %in% dexseq.res.significantID,]
  HELPER_SAVE_DATA_FRAME(as.data.frame(resDataFrame),paste0(outputDir,"/DEA/RES_SIGNIFICANT_GENES.csv") )
  
  dexseq.res.significant.glycoID = unique(dexseq.res.significant$groupID[dexseq.res.significant$groupID %in% geneList])
  dexseq.res.significant.glyco = dexseq.res.significant[dexseq.res.significant$groupID %in% dexseq.res.significant.glycoID,]
  HELPER_SAVE_DATA_FRAME(as.data.frame(dexseq.res.significant.glyco),paste0(outputDir,"/DEA/RES_SIGNIFICANT_GLYCO_GENES.csv") )
  return (list(significantID.all = dexseq.res.significantID,significantID.glyco = dexseq.res.significant.glycoID))
}

HELPER_SAVE_ORCA = function (graphics, outputDir, fileName){
  current = getwd()
  setwd(outputDir)
  orca(graphics,fileName)
  setwd(current)
}