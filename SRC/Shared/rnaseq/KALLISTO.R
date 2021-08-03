#esto no va mas--> lo movi a htseq.r
#*******************************************************************************************************
#Kallisto + DESEQ2 Gene level expression
#*******************************************************************************************************
outputDir.kallisto =  paste0(project.resultsPath, "KALLISTO/")
CreateGroupFolder (outputDir =outputDir.kallisto,project.groupsToCompare)

files = paste0(project.resourcePath,"kallisto/",samplesList$FILE , "/abundance.h5")

kallisto.geneLevel.dds = DESeq2_tx_import(filesPath = files,transcriptAnnotation ,type = "kallisto",  geneLevel= TRUE, conditionDataFrame = samplesList, formula = ~GROUP, coresNumber = 4)
kallisto.geneLevel.vsd = DESeq2_VST (kallisto.geneLevel.dds)

kallisto.rawCountMatrix = assay(kallisto.geneLevel.dds)
colnames(kallisto.rawCountMatrix) = samplesList$FILE


kallisto.rawCountMatrix.annotation =  HELPER_RowNamesAsFirstColumn(kallisto.rawCountMatrix, "ENSENMBL_GENE_ID")
kallisto.rawCountMatrix.annotation = ANN_MergeAnnotationToDataFrame(kallisto.rawCountMatrix.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")

HELPER_SAVE_DATA_FRAME(kallisto.rawCountMatrix.annotation, paste0(project.resultsPath, "KALLISTO/CountMatrix/rawCountMatrixAnnotated.csv"))

DA_RNASEQ_ALL(rawCountMatrix = kallisto.rawCountMatrix,
       deseq.dds = kallisto.geneLevel.dds,
       dseq.vsd = kallisto.geneLevel.vsd,
       samplesNames = samplesList$SAMPLE,
       samplesGroups = samplesList$GROUP,
       tsne.perplexity = 2,
       KMeans.kmax = length(levels(samplesList$GROUP))*3,
       KMeans.optimalNumberOfClusters =  length(levels(samplesList$GROUP)),
       KMeans.maxIteration = settings.kmeans.numberOfinteractions,
       KMeans.numbersOfRandomCenters = settings.kmeans.numbersOfRandomCenters,
       groups.colorArray = settings.color,
       outputDir = outputDir.kallisto)


#compare groups
for (i in 1:length(project.groupsToCompare)){
  contrast =project.groupsToCompare[[i]]

  outputDirGroup = paste0(outputDir.kallisto,"DEA/",project.groupsToCompare[[i]][2],"_VS_",project.groupsToCompare[[i]][3],"/")
  
  #NOTA aca la hipotesis nula es que log2FC <=1 y la alternativa >1
  kallisto.geneLevel.res = DESeq2_ExtractResults (kallisto.geneLevel.dds,filterFun = ihw, coresNumber = 4, pAdjustMethod = "BH", adjustedPValue = project.settings.DEA.padj, contrast = contrast,  format ="DataFrame",lfcThreshold = project.settings.DEA.H0.Log2FC)
  kallisto.geneLevel.shrunken <- lfcShrink(kallisto.geneLevel.dds,coef = resultsNames(kallisto.geneLevel.dds)[2], lfcThreshold = 0)  #SOLO para plotear y raankear, despues da lo mismo

  kallisto.geneLevel.res.annotation = HELPER_RowNamesAsFirstColumn(kallisto.geneLevel.res, "ENSENMBL_GENE_ID")
  kallisto.geneLevel.res.annotation = ANN_MergeAnnotationToDataFrame(kallisto.geneLevel.res.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")

  kallisto.geneLevel.shrunken.annotation = HELPER_RowNamesAsFirstColumn(kallisto.geneLevel.shrunken, "ENSENMBL_GENE_ID")
  kallisto.geneLevel.shrunken.annotation = ANN_MergeAnnotationToDataFrame(kallisto.geneLevel.shrunken.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")

  
  DESeq2_Result_ALL_Stats(DESeqResult = kallisto.geneLevel.res, 
                          DESeqResult.shrunken =  kallisto.geneLevel.shrunken, 
                          DESeqResult.annotation = kallisto.geneLevel.res.annotation,
                          padj = project.settings.DEA.padj,
                          log2FCCutOff = project.settings.DEA.cutoff.log2FC, 
                          outputDir = outputDirGroup)
  
  HELPER_FILTER_SAVE_PFI (DESeqResult.annotation = kallisto.geneLevel.res.annotation,
                               geneList = geneFilteringList,
                               padj = project.settings.DEA.padj, 
                               log2FC = project.settings.DEA.cutoff.log2FC, 
                               outputDir = outputDirGroup 
  )

  

  #*******************************************************************************************************
  #GO   
  #*******************************************************************************************************
  #goseq
  #Nota aca como pusimos lfcthreshold =1, el shrunken noda padj, entonces pongo el res normal, pero no se si esta bien
  ONTHO_GO_MOUSE_ALL (kallisto.geneLevel.shrunken,n = 20,genome = "mm10", geneID = "ensGene", padj = project.settings.DEA.padj,log2FC_UP = project.settings.DEA.cutoff.log2FC, log2FC_DOWN = -project.settings.DEA.cutoff.log2FC, categories = c("GO:BP"), outputDir  = paste0(outputDirGroup,"GO/"))
  
  
  #*******************************************************************************************************
  #GSEA  
  #*******************************************************************************************************
  #fgsea
  ranks = ONTO_GSEA_RANK_GENES(deseq.res.annotation = kallisto.geneLevel.shrunken.annotation, 
                               rankType ="log2FC",
                               outputDir =  paste0(outputDirGroup,"GSEA/"))
  dbs = ONGTHO_GSEA_GET_MOUSE_DB()
  for (i in 1:length(dbs)){
    fgseaRes = ONTHO_GSEA_MOUSE_Calculate(ranks, misegdbType = dbs[i], outputDir = paste0(outputDirGroup,"GSEA/"), coresNumber = setttings.coresNumber)
    
    if (table(fgseaRes$result$padj < project.settings.DEA.padj)[2]<50){#dibujamos solos que tienen pocos 
      ONTHO_GSEA_Plot_Result(fgseaRes = fgseaRes$result, padj = project.settings.DEA.padj, title =  dbs[i], outputDir = paste0(outputDirGroup,"GSEA/"), dbName = dbs[i])
    }
    
  }
  #*******************************************************************************************************
  #KEGG 
  #*******************************************************************************************************
  search_kegg_organism('mmu', by='kegg_code')
  sigGenes <- kallisto.geneLevel.shrunken.annotation$entrezgene[ kallisto.geneLevel.shrunken.annotation$padj < project.settings.DEA.padj & 
                                                       !is.na(kallisto.geneLevel.shrunken.annotation$padj) &
                                                       abs(kallisto.geneLevel.shrunken.annotation$log2FoldChange) > project.settings.DEA.cutoff.log2FC ]
  sigGenes <- na.exclude(sigGenes)
  kk <- enrichKEGG(gene = sigGenes, organism = 'mmu')
  HELPER_SAVE_DATA_FRAME(as.data.frame(kk),paste0(outputDirGroup,"KEGG/keggPathwaysALL.csv"))
  
  kegg.significant = as.data.frame(kk)

  kegg.significant =  kegg.significant$ID[kegg.significant$p.adjust<project.settings.DEA.padj]
  
  
  kegg.log2fc <- kallisto.geneLevel.res.annotation$log2FoldChange  
  names(kegg.log2fc) <- kallisto.geneLevel.res.annotation$entrezgene
  tempFolder = getwd()
  setwd(paste0(outputDirGroup,"KEGG/"))
  for (i in 1:length(kegg.significant)){
    try(
        pathview(gene.data = kegg.log2fc, 
                 pathway.id = kegg.significant[i], 
                 species = "mmu", 
                 limit = list(gene=5, cpd=1),
                 out.suffix = project.name
        )
    )       
  }
  setwd(tempFolder)
  
}
