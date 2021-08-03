#Htseq y kallisto Tienen que estar todo en 1, lo mismo para TCGA
#TODO tcgabiolinks o rtcgao arch4 o xena como traigo la data y clinical firehouse
#agregar para kalisto isoform package 
#TODO de tcga hay que homogeneizar la sample list para poder tener un solo archivo de RNASEQ
#************************************************************
#TODO GO y GSEA KEGG independiente de la especie 
#ver el merge para que quede siempre el symbol primero y quede del mismo tamaño que la matrix ANN_MergeAnnotationToDataFrame(geneAnnotation.immune,htseqFPKM.rawCountMatrix.annotation,"ensembl_gene_id","ENSENMBL_GENE_ID", keepAll =  "Rigth")
#************************************************************
#TODO formula con varios factores y results
#TODO como ploteo con varios factores hacer varia pruebas con datos completos, con missing values y con control sin otros datos
#TODO cambiar en todo para guardar tsv
#TODO El gene ratio tienen solo para ratones y hay que poner en humanos u otro
#*******************************************************************************************************
#GENE LEVEL
#*******************************************************************************************************
# HTSeq + DESeq2  +++   KALISTO --> hay que sacar el nombre de htseq a las variables
#*******************************************************************************************************
if (project.rnaseq.runHTSEQ){
  outputDir.htseq =  paste0(project.resultsPath, "HTSEQ/")
  
  htseq.rawCountMatrix = DESeq2_MergeHTSeqExpressionFiles(paste0(project.resourcePath,"htseq"),samplesList$HTSEQ_FILE)
  htseq.dds = DESeq2_Matrix_Import(htseq.rawCountMatrix, samplesList, project.formula, coresNumber = setttings.coresNumber)
  
}else{
  outputDir.htseq =  paste0(project.resultsPath, "KALLISTO/")
  
  
  files = paste0(project.resourcePath,"kallisto/",samplesList$FILE , "/abundance.h5")
  
  htseq.dds = DESeq2_tx_import(filesPath = files,transcriptAnnotation ,type = "kallisto",  geneLevel= TRUE, conditionDataFrame = samplesList, formula =project.formula, coresNumber = setttings.coresNumber)
  
  
  htseq.rawCountMatrix = assay(htseq.dds)
  colnames(htseq.rawCountMatrix) = samplesList$FILE
  
}
CreateGroupFolder(outputDir =outputDir.htseq,project.groupsToCompare)

htseq.rawCountMatrix.annotation =  HELPER_RowNamesAsFirstColumn(htseq.rawCountMatrix, "ENSENMBL_GENE_ID")
htseq.rawCountMatrix.annotation = ANN_MergeAnnotationToDataFrame(htseq.rawCountMatrix.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")
HELPER_SAVE_DATA_FRAME(htseq.rawCountMatrix.annotation, paste0(outputDir.htseq, "CountMatrix/rawCountMatrixAnnotated.csv"))

htseq.vsd = DESeq2_VST(htseq.dds)
HELPER_SAVE_DATA_FRAME(as.data.frame(assay(htseq.vsd)), paste0(outputDir.htseq, "CountMatrix/vsd.csv"))

if (project.rnaseq.run.da){
    DA_RNASEQ_ALL(rawCountMatrix = htseq.rawCountMatrix,
           deseq.dds = htseq.dds,
           dseq.vsd = htseq.vsd,
           samplesNames = project.sample.name,
           samplesGroups = project.sample.group,
           tsne.perplexity = 2,
           KMeans.kmax =  length(levels(project.sample.group))*3,
           KMeans.optimalNumberOfClusters = length(levels(project.sample.group)),
           KMeans.maxIteration = settings.kmeans.numberOfinteractions,
           KMeans.numbersOfRandomCenters = settings.kmeans.numbersOfRandomCenters,
           groups.colorArray = settings.color,
           outputDir = outputDir.htseq)
    
}
  #compare groups
for (i in 1:length(project.groupsToCompare)){
    contrast =project.groupsToCompare[[i]]
    for (j in 1:length(project.settings.DEA.H0.Log2FC)){
        tempLog2FC = project.settings.DEA.H0.Log2FC[j]
        outputDirGroup = paste0(outputDir.htseq,"DEA/",project.groupsToCompare[[i]][2],"_VS_",project.groupsToCompare[[i]][3],"_Log2FC_",tempLog2FC,"/")
       
        
        htseq.res = DESeq2_ExtractResults (htseq.dds,filterFun = ihw, coresNumber = setttings.coresNumber, pAdjustMethod = "BH", adjustedPValue = project.settings.DEA.padj,contrast,  format ="DataFrame", lfcThreshold = tempLog2FC)
        htseq.shrunken <- lfcShrink(htseq.dds,coef = resultsNames(htseq.dds)[2], lfcThreshold = tempLog2FC)  #SOLO para plotear y raankear, despues da lo mismo
        
        htseq.res.annotation = HELPER_RowNamesAsFirstColumn(htseq.res, "ENSENMBL_GENE_ID")
        htseq.res.annotation = ANN_MergeAnnotationToDataFrame(htseq.res.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")
        
        hstseq.shrunken.annotation = HELPER_RowNamesAsFirstColumn(htseq.shrunken, "ENSENMBL_GENE_ID")
        hstseq.shrunken.annotation = ANN_MergeAnnotationToDataFrame(hstseq.shrunken.annotation,geneAnnotation,"ENSENMBL_GENE_ID","ensembl_gene_id", keepAll =  "Left")
        
        
        HELPER_FILTER_SAVE_PFI (DESeqResult.annotation = htseq.res.annotation,
                                geneList = geneFilteringList,
                                padj = project.settings.DEA.padj, 
                                log2FC = project.settings.DEA.cutoff.log2FC, 
                                outputDir = outputDirGroup 
        )
        
        DESeq2_Result_ALL_Stats(DESeqResult = htseq.res, 
                                DESeqResult.shrunken =  htseq.shrunken, 
                                DESeqResult.annotation = htseq.res.annotation,
                                padj = project.settings.DEA.padj,
                                log2FCCutOff = project.settings.DEA.cutoff.log2FC, 
                                outputDir = outputDirGroup
        )
      if (project.rnaseq.run.go){    
          if (project.specie =="MOUSE"){
            dir.create(file.path(paste0(outputDirGroup,"GO")), recursive = TRUE, showWarnings = FALSE)
            dir.create(file.path(paste0(outputDirGroup,"GSEA")), recursive = TRUE, showWarnings = FALSE)
            dir.create(file.path(paste0(outputDirGroup,"KEGG")), recursive = TRUE, showWarnings = FALSE)
             #*******************************************************************************************************
              #GO   
              #*******************************************************************************************************
              #goseq
              ONTHO_GO_MOUSE_ALL (htseq.shrunken,n = 20,genome = "mm10", geneID = "ensGene", padj = project.settings.DEA.padj,log2FC_UP = project.settings.DEA.cutoff.log2FC, log2FC_DOWN = -project.settings.DEA.cutoff.log2FC, categories = c("GO:BP"), outputDir  = paste0(outputDirGroup,"GO/"))
          
          
              #*******************************************************************************************************
              #GSEA  
              #*******************************************************************************************************
              #fgsea
              ranks = ONTO_GSEA_RANK_GENES(deseq.res.annotation = hstseq.shrunken.annotation,
                                           rankType ="log2FC",
                                           outputDir =  paste0(outputDirGroup,"GSEA/"))
              dbs = ONGTHO_GSEA_GET_MOUSE_DB()
              for (i in 1:length(dbs)){
                fgseaRes = ONTHO_GSEA_MOUSE_Calculate(ranks, misegdbType = dbs[i], outputDir = paste0(outputDirGroup,"GSEA/"), coresNumber = setttings.coresNumber)
          
                if (table(fgseaRes$result$padj<project.settings.DEA.padj)[2]<50){#dibujamos solos que tienen pocos 
                  ONTHO_GSEA_Plot_Result(fgseaRes = fgseaRes$result, padj = project.settings.DEA.padj, title =  dbs[i], outputDir = paste0(outputDirGroup,"GSEA/"), dbName = dbs[i])
                }
                
              }
              #*******************************************************************************************************
              #KEGG 
              #*******************************************************************************************************
              search_kegg_organism('mmu', by='kegg_code')
              sigGenes <- hstseq.shrunken.annotation$entrezgene[ hstseq.shrunken.annotation$padj < project.settings.DEA.padj & 
                                                                   !is.na(hstseq.shrunken.annotation$padj) &
                                                                   abs(hstseq.shrunken.annotation$log2FoldChange) > project.settings.DEA.cutoff.log2FC ]
              sigGenes <- na.exclude(sigGenes)
              kk <- enrichKEGG(gene = sigGenes, organism = 'mmu')
              HELPER_SAVE_DATA_FRAME(as.data.frame(kk),paste0(outputDirGroup,"KEGG/keggPathwaysALL.csv"))
              
              kegg.significant = as.data.frame(kk)
              kegg.significant.description = kegg.significant$Description[kegg.significant$p.adjust < project.settings.DEA.padj]
              kegg.significant = kegg.significant$ID[kegg.significant$p.adjust < project.settings.DEA.padj]
          
              
              kegg.log2fc <- htseq.res.annotation$log2FoldChange  
              names(kegg.log2fc) <- htseq.res.annotation$entrezgene
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
      }
    }
}
















