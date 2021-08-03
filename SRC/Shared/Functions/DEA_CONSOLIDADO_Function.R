
collapseGenesDataFrame = function (resDataFrame){
  resDataFrame = resDataFrame[order(resDataFrame$external_gene_name),]
  
  #as numeric
  matrix = resDataFrame[,4:dim(resDataFrame)[2]]
  for (i in 1:dim(matrix)[2]){
    matrix[,i] = as.numeric(matrix[,i] ) 
  }
  #collapse duplicated genes
  collapsedMatrix = matrix
  collapsedMatrix$geneName = resDataFrame$external_gene_name
  collapsedMatrix= collapsedMatrix %>% 
    group_by(geneName) %>% 
    summarise_all( mean, na.rm = TRUE)
  dim(collapsedMatrix)
  
  #get the families for each gene
  geneFamily = resDataFrame[collapsedMatrix$geneName %in% resDataFrame$external_gene_name,]
  geneFamily = geneFamily[!duplicated(geneFamily$external_gene_name),]
  dim(geneFamily)
  geneFamily = geneFamily[order(match(geneFamily$external_gene_name,collapsedMatrix$geneName)),]
  collapsedMatrix$family = geneFamily$Group
  
  return (collapsedMatrix)
}
#ENSEMBL_GENE_ID
readAndMergeFiles = function(resultPath, geneList, outputDir, filePrefix = "HS", mergeByGeneName = FALSE){
  
  #get all tcga significat gylco genes from tcga Result
  tcgaFiles = dir(resultPath, recursive = TRUE)
  tcgaFiles = tcgaFiles[grep("*RES_SIGNIFICANT_GLYCO_GENES.csv", tcgaFiles)]
  removeDexseq = tcgaFiles[grep("*DEXSeq",tcgaFiles)]
  tcgaFiles = tcgaFiles[!tcgaFiles %in%removeDexseq]
  
  tcgaFiles.logfc1 = tcgaFiles[grep(".*_Log2FC_1.*", tcgaFiles)]
  names.logfc1 = gsub("/DEA/RES_SIGNIFICANT_GLYCO_GENES.csv", "",tcgaFiles.logfc1)
  names.logfc1 = gsub("HTSEQ/DEA/", "",names.logfc1)
  names.logfc1 = gsub("^COAD_.*/", "",names.logfc1)
  names.logfc1  = gsub("_Log2FC_1", "",names.logfc1)
  names.logfc1
  
  tcgaFiles.logfc0 = tcgaFiles[!tcgaFiles %in%tcgaFiles.logfc1]
  names.logfc0 = gsub("/DEA/RES_SIGNIFICANT_GLYCO_GENES.csv", "",tcgaFiles.logfc0)
  names.logfc0 = gsub("HTSEQ/DEA/", "",names.logfc0)
  names.logfc0 = gsub("^COAD_.*/", "",names.logfc0)
  names.logfc0
  
  tcgaRes.logfc0 = geneList
  
  #merge files
  for (i in 1:length(tcgaFiles.logfc0)){
    a = HELPER_LoadTSV(paste0(resultPath,tcgaFiles.logfc0[i]))

    if (mergeByGeneName){
      a = data.frame(external_gene_name = a$external_gene_name,log2FoldChange= a$log2FoldChange)
      colnames(a) = c("external_gene_name", names.logfc0[i])
      if (is.na(tcgaRes.logfc0)){
        tcgaRes.logfc0 = a
      }else{
        tcgaRes.logfc0 = merge(tcgaRes.logfc0, a, by.x = "external_gene_name", by.y = "external_gene_name", all.x = TRUE)
      }
    }else{
      a = data.frame(ensembl_gene_id = a$ensembl_gene_id,log2FoldChange= a$log2FoldChange)
      colnames(a) = c("ensembl_gene_id", names.logfc0[i])
      if (is.na(tcgaRes.logfc0)){
        tcgaRes.logfc0 = a
      }else{
        tcgaRes.logfc0 = merge(tcgaRes.logfc0, a, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
      }
    }
  }
  #View(res.logfc0)
  HELPER_SAVE_DATA_FRAME(tcgaRes.logfc0, paste0(project.resultsPath,filePrefix,"_glyco.csv"))
  
  tcgaRes.logfc1 = geneListHS
  for (i in 1:length(tcgaFiles.logfc1)){
    a = HELPER_LoadTSV(paste0(resultPath,tcgaFiles.logfc1[i]))
    if (mergeByGeneName){
      a = data.frame(external_gene_name = a$external_gene_name,log2FoldChange= a$log2FoldChange)
      colnames(a) = c("external_gene_name", names.logfc0[i])
      if (is.na(tcgaRes.logfc1)){
        tcgaRes.logfc1 = a
      }else{
        tcgaRes.logfc1 = merge(tcgaRes.logfc0, a, by.x = "external_gene_name", by.y = "external_gene_name", all.x = TRUE)
      }
    }else{
      a = data.frame(ensembl_gene_id = a$ensembl_gene_id,log2FoldChange= a$log2FoldChange)
      colnames(a) = c("ensembl_gene_id", names.logfc1[i])
      if (is.na(tcgaRes.logfc1)){
        tcgaRes.logfc1 = a
      }else{
        tcgaRes.logfc1 = merge(tcgaRes.logfc1, a, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
      }
    }
  }
  #View(res.logfc1)
  HELPER_SAVE_DATA_FRAME(tcgaRes.logfc1, paste0(project.resultsPath,filePrefix,"glycoLog2FC1.csv"))
  
  return (list(H0log2FC0 = tcgaRes.logfc0,H0log2FC1 = tcgaRes.logfc1))
}
