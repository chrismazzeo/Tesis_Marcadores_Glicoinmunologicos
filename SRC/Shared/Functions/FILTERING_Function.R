#********************
# Filtering  
#*******************
FILTERING_DESEQ2_RES_ByAdjPvalue = function (res, padj){
  temp = res[!is.na(res$padj) & res$padj < padj,]
  
  return (temp)
  
}

FILTERING_ByGeneList = function (sourceDataFrame,sourceColumn,filterColumn){
  temp = sourceDataFrame[sourceColumn %in% filterColumn,]
  return (temp)
}

#TODO esto hay que mejorarlo para que filtre cualquier campo que le pase
#Filter Sample Type & remove duplicated  & sort
FILTERING_Samples = function (inputFile, outputPath = NULL){
  
  samples = loadTSV(inputFile)
  #seleccionamos que tipo de tejido
  samples = samples [samples$`Sample Type` == "Primary Tumor" | samples$`Sample Type` == "Solid Tissue Normal",]
  #removemos top down
  samples = samples %>% filter(str_detect(samples$`Sample ID`,"-01A") | str_detect(samples$`Sample ID`,"-11A"))
  
  #removemos replicas
  samples = samples[!duplicated(samples$`Sample ID`),]
  
  #ordenamos
  samples = samples[order(samples$`Sample Type`, samples$`Sample ID`),]
  
  
  saveTSV(samples,outputPath)
  
  
  return (samples)
}
