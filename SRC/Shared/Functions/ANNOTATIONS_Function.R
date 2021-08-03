

list.of.packages <- c("readr","biomaRt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocParallel")
#biocLite("apeglm")
#biocLite("ReportingTools")
#biocLite("org.Hs.eg.db")

library(readr)
library(biomaRt)

#***************************************************************
#Annotations extra
#***************************************************************
ANN_GetEnsemblSpecies = function(){
  mart = useEnsembl(biomart="ensembl")
  useEnsembl('ENSEMBL_MART_ENSEMBL')
  return (listDatasets(mart))
}

#****
# return all typs of ids which can be retrieve for a specie
#****
ANN_GetEnsemblAttributeList = function (specie){
  mart = useEnsembl(biomart="ensembl", dataset=specie)
  return (listAttributes(mart))
}


#***********************
# ensemblIDs = genes IDs to retrieve annotation infor
# ensemblIDType  = ensembl gene id, 
# annotationIDList = what type of attributes you want to retrieve

#annotationIDList = c('hgnc_symbol','gene_biotype','superfamily','family','family_description','phenotype_description','chromosome_name','start_position','end_position','description','strand','band')
#geneOrthology = c('mmusculus_homolog_ensembl_gene','mmusculus_homolog_associated_gene_name','mmusculus_homolog_orthology_type')
#***********************
ANN_GetAnnotationAll = function(specie = NULL, ensemblIDType = "ensembl_gene_id_version",attributes,collapseInRow = TRUE, uniqueRows = TRUE){
  
  if (is.null(specie)){
    print ("Falta la especie")
    return()
  }
  mart = useEnsembl(biomart="ensembl", dataset=specie)

  annotations <- getBM(attributes= attributes,
                       mart= mart,
                       uniqueRows = uniqueRows)
  if (collapseInRow){
    # annotations = ANN_CollapseAnnotationsInOneRow(annotations)  
    annotations = annotations[!duplicated(annotations$ensembl_gene_id_version),]
  }
  
  
  return(annotations)
}

ANN_GetAnnotation = function(specie = "hsapiens_gene_ensembl", ensemblIDType = "ensembl_gene_id_version",filterIDs, attributes, uniqueRows = TRUE){
  
  mart = useEnsembl(biomart="ensembl", dataset=specie)
  
  
  annotations <- getBM(filters= ensemblIDType,
                         values = filterIDs,
                         attributes= attributes,
                         mart= mart,
                         uniqueRows = uniqueRows)
  
  
  return(annotations)
}

#*******
# This will collapse all the results in one row per gene ID
#******
ANN_CollapseAnnotationsInOneRow = function (annotations){
  #annotations = order(annotations, annotations[1], decreasing = FALSE)
  annotations = annotations[order(annotations[1], decreasing  = FALSE),]
 
    mergedDataFrame = annotations[1,]
    colnames(mergedDataFrame) = colnames(annotations)
    j = 1
    for (i in 2:dim(annotations)[1]){
      if (annotations[i,1]  == mergedDataFrame[j,1]){ #merge
        for (k in 2:dim(annotations)[2]){ #veo que columna tengo que mergear
          if  (!is.na(annotations[i,k]) & (annotations[i,k] != mergedDataFrame[j,k]) ){  #son distintos
            if (!is.na(annotations[i,k]) & annotations[i,k] != ""){
              mergedDataFrame[j,k] = paste0(mergedDataFrame[j,k], " - ", annotations[i,k])
            }
          }
          else{
            mergedDataFrame[j,k] = annotations[i,k]
          }
        }
      }
      else{ #es nuevo y lo agrego
        mergedDataFrame = rbind(mergedDataFrame,annotations[i,])
        j = j+1
      }
    }
    #annotations = mergedDataFrame
  return (mergedDataFrame)
}


ANN_MergeAnnotationToDataFrame = function (resultDataFrame, annotationsDataFrame,leftID,rigthID, keepAll ="none"){
 
  
  if (keepAll == "None"){
    tempRes= merge(x=resultDataFrame, y=annotationsDataFrame,  by.x = leftID, by.y =rigthID )
    return (tempRes)
  }
  if (keepAll == "Left"){
    tempRes= merge(x=resultDataFrame, y=annotationsDataFrame,  by.x = leftID, by.y =rigthID, all.x = TRUE )
    return (tempRes)
  }
  
    if (keepAll == "Rigth"){
      tempRes= merge(x=resultDataFrame, y=annotationsDataFrame,  by.x = leftID, by.y =rigthID, all.y = TRUE )
      return (tempRes)
  }
  #countMatrix = cbind(ensembl_gene_id = countMatrix$ensembl_gene_id, symbol = countMatrix$hgnc_symbol,countMatrix[,2:(dim(countMatrix)[2]-1)])
 
}

ANN_RemoveStringAfterDot = function(data){
  strings = gsub("\\..*","",data)  
  return (strings)
}

