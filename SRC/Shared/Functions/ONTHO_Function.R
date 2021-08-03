#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#go-enrichment-analysis

#gene onthology
#solo vemos que GO estan enriquecidos (o sea over)
ONTHO_GO_ORGANISM =function (){
  return (supportedOrganisms())
}
ONTHO_GO_GENE_LIST = function (deseq.res, padj, log2FC_UP, log2FC_DOWN){
  #de la lista de DE hay que sacar los que tengan padj = NA, de ahi armar 2 columnas 1 con todos los genes y otra con con los que significativos 1 significativo, 0 no significativo
  universe = rownames(deseq.res)[!is.na(deseq.res$padj)]
  myInterestedGenes <- rownames(deseq.res)[!is.na(deseq.res$padj) & deseq.res$padj<padj & (deseq.res$log2FoldChange> log2FC_UP | deseq.res$log2FoldChange< log2FC_DOWN)]
  length(myInterestedGenes)
  geneList <- factor(as.integer(universe %in% myInterestedGenes))
  names(geneList) <- universe
  str(geneList)
  return (geneList)
}

ONTHO_GOSEQ_RESULT = function(geneList, genome, geneID, padj = 0.05, outputDir = NULL,categories = c("GO:BP")){
  #GOseq is a method to conduct Gene Ontology (GO) analysis suitable for RNA-seq data as it accounts for the gene length bias in detection of over-representation
  #******************************
  #goseq OK
  #******************************
  #aca usamos como lista de genes solo los que son significativos
  #DEGenes
  #bias.data giving the numeric value of the DE bias  being accounted for (usually the gene length or number of counts)
  #pwf giving the genes value on the probability weighting function.
  pdf(paste0(outputDir,"fit.pdf"))
  pwf=nullp(geneList,genome,geneID) #devuelve 3 columnas, 
  dev.off()
  
  
  
  #Onthology (es la categoria a la que pertenece cada categoria GO) BP biological process, CC celular components, MP molecular functions
  #ojo que el pvalue aca no tienen correccion por multiples correcciones, se hace abajo esto
  goResults = goseq(pwf,genome,geneID,test.cats=categories)

  #interpretacion, corregimos el pvalue con p.adjust() para obtener un padj
  enriched.GO.up = data.frame(goResults, padj = p.adjust(goResults$over_represented_pvalue, method="BH"))
  enriched.GO.up = enriched.GO.up[enriched.GO.up$padj <padj,]
  
  goseq.result.over = data.frame(GOID =enriched.GO.up$category,
                              TERM = enriched.GO.up$term,
                              DEFINITION =  Definition(enriched.GO.up$category),
                              Ontology =  enriched.GO.up$ontology,
                              numDEInCat = enriched.GO.up$numDEInCat,
                              numInCat = enriched.GO.up$numInCat,
                              over_represented_pvalue = enriched.GO.up$over_represented_pvalue,
                              padj = enriched.GO.up$padj
  )

  HELPER_SAVE_DATA_FRAME(goseq.result.over,paste0(outputDir,"EnrichedOver.csv"))
 
   
  return (goseq.result.over)
}
ONTHO_GO_PLOT = function (goseq.result.over, n = 10, outputDir){
  graphicsOver =  goseq.result.over %>% 
    top_n(n = n, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc, 
               y=TERM, 
               colour=over_represented_pvalue, 
               size=numDEInCat)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
  

  
 print(graphicsOver)

 HELPER_SAVE_PDF(graphicsOver,paste0(outputDir,"over.pdf") )
}
ONTHO_GO_MOUSE_ALL = function (deseq.shurken,n = 20 ,genome = "mm10", geneID = "ensGene", padj = 0.05,log2FC_UP = 1, log2FC_DOWN = -1, categories = c("GO:CC","GO:BP","GO:MF"), outputDir ){

  geneList = ONTHO_GO_GENE_LIST(deseq.shurken,padj = padj, log2FC_UP = log2FC_UP, log2FC_DOWN = log2FC_DOWN)
  go.res = ONTHO_GOSEQ_RESULT(geneList, genome = genome, geneID = geneID, padj = padj, outputDir = outputDir, categories )
  ONTHO_GO_PLOT(go.res, n = n, outputDir = outputDir)
  
}

ONTHO_GO_MOUSE_ALL = function (deseq.shurken,n = 20 ,genome = "mm10", geneID = "ensGene", padj = 0.05,log2FC_UP = 1, log2FC_DOWN = -1, categories = c("GO:CC","GO:BP","GO:MF"), outputDir ){
  
  
  geneList = ONTHO_GO_GENE_LIST(deseq.shurken,padj = padj, log2FC_UP = log2FC_UP, log2FC_DOWN = log2FC_DOWN)
  go.res = ONTHO_GOSEQ_RESULT(geneList, genome = genome, geneID = geneID, padj = padj, outputDir = outputDir, categories )
  ONTHO_GO_PLOT(go.res, n = n, outputDir = outputDir)
  
}

#GSE
#aca usamos todos los genes pero los rankeamos
#sacamos los que no tengan entreID
#si usamos el shrunken podemos rankear por log2fc o -log10(pvalue)* sign(log2FC)
#si no usamos el shrunken, rankeamos por stats

#  The analysis is performed by:
#ranking all genes in the data set
#identifying the rank positions of all members of the gene set in the ranked data set
#calculating an enrichment score (ES) that represents the difference between the observed rankings and that which would be expected assuming a random rank distribution.
ONTO_GSEA_RANK_GENES = function(deseq.res.annotation, rankType = "log2FC",outputDir = NULL){
  fgsea.data = deseq.res.annotation[!is.na(deseq.res.annotation$entrezgene),]
  
  if (rankType == "log2FC"){ #shrunken
    if("stats" %in% colnames(deseq.res.annotation))
    {
      print ("hay que pasar el shrunken")
    }
    else{
      ranks  = fgsea.data$log2FoldChange
  
    }
  }
  else if (rankType == "stats"){
    ranks  = fgsea.data$stats
  }
  else if (rankType == "pvalue"){#shrunken
    if("stats" %in% colnames(deseq.res.annotation))
    {
      print ("hay que pasar el shrunken")
    }
    else{
      ranks  = -log10(pvalue) * sign(fgsea.data$log2FoldChange)
    }
  }
  
  names(ranks) = fgsea.data$entrezgene
  ranks = sort(ranks, decreasing = T) 
  head(ranks)
  if (!is.null(outputDir)){
    pdf(paste0(outputDir,"rankby_",rankType,".pdf"))
    barplot(ranks)
    dev.off()
  }
  barplot(ranks)
  return (ranks)
}
ONGTHO_GSEA_GET_MOUSE_DB = function (){
    return (list("Human","Curated","Motif","Computational","GO","Oncogenic","Immunologic"))
}

ONTHO_GSEA_MOUSE_Calculate = function (ranks, misegdbType = "Human",outputDir = NULL, coresNumber = 1){
  #para mouse tenemos que usar los homologos por eso cargamos una base de datos homologoa de humanos
  load ("./ExternalLibs/wehi/mouse_H_v5p2.rdata") # es el equivalente de mouse de misigDB // http://bioinf.wehi.edu.au/software/MSigDB/index.html
  #aca hay distintos patwsway que se pueden usar y depende del archivo que bajemos
  

  
  if (misegdbType == "Human"){
    load ("./ExternalLibs/wehi/mouse_H_v5p2.rdata")
    pathways = Mm.H  
  }

  else if (misegdbType == "Curated"){
    load ("./ExternalLibs/wehi/mouse_c2_v5p2.rdata")
    pathways = Mm.c2  
  }
  else if (misegdbType == "Motif"){
    load ("./ExternalLibs/wehi/mouse_c3_v5p2.rdata")
    pathways = Mm.c3 
  }
  else if (misegdbType == "Computational"){
    load ("./ExternalLibs/wehi/mouse_c4_v5p2.rdata")
    pathways = Mm.c4 
  }
  else if (misegdbType == "GO"){
    load ("./ExternalLibs/wehi/mouse_c5_v5p2.rdata")
    pathways = Mm.c5 
  }
  else if (misegdbType == "Oncogenic"){
    load ("./ExternalLibs/wehi/mouse_c6_v5p2.rdata")
    pathways = Mm.c6  
  }
  else if (misegdbType == "Immunologic"){
    load ("./ExternalLibs/wehi/mouse_c7_v5p2.rdata")
    pathways = Mm.c7  
  }

  names(pathways)

  if (coresNumber >1){
    BPPARAM = MulticoreParam(workers=coresNumber)
  }
  #ES – enrichment score, same as in Broad GSEA implementation  ES >0 es enriquecido y ES <0 esta deseqnrquicido
  #NES – enrichment score normalized to mean enrichment of random samples of the same size;
  fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize=15, maxSize = 500, nperm=1000,BPPARAM = BPPARAM) 
  #The warning produced indicates that there are few genes that have the same fold change and so are ranked equally. fgsea with arbitrarily order determine which comes first in the ranked list. 
  #0 -1% aprox   - As long as this number is small it shouldn’t significantly effect the results. If the number is large something is suspicious about the fold change results.
  
  fgseaRes = fgseaRes[order(padj, -abs(NES)),]
  fgseaRes.save = fgseaRes

  
  for (i in 1: length(fgseaRes.save$leadingEdge)){
    collapesed = as.character(paste( unlist(fgseaRes.save$leadingEdge[i]), collapse=' - '))  
    fgseaRes.save$ENTREZ_GENES_IDS_IN[i] = collapesed
    
    
  }
  fgseaRes.save = fgseaRes.save[,-8]

  
  HELPER_SAVE_DATA_FRAME(as.data.frame(fgseaRes.save), paste0(outputDir,misegdbType,".csv"), rowNames = FALSE)

  return (list (result = fgseaRes, pathways = pathways))
}

ONTHO_GSEA_HUMAN_Calculate = function (){
  #para humanos #https://stephenturner.github.io/deseq-to-fgsea/#using_the_fgsea_package
  
  if (coresNumber >1){
    BPPARAM = MulticoreParam(workers=coresNumber)
  }
  #ES – enrichment score, same as in Broad GSEA implementation  ES >0 es enriquecido y ES <0 esta deseqnrquicido
  #NES – enrichment score normalized to mean enrichment of random samples of the same size;
  fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize=15, maxSize = 500, nperm=1000,BPPARAM = BPPARAM) 
  #The warning produced indicates that there are few genes that have the same fold change and so are ranked equally. fgsea with arbitrarily order determine which comes first in the ranked list. 
  #0 -1% aprox   - As long as this number is small it shouldn’t significantly effect the results. If the number is large something is suspicious about the fold change results.
  
  fgseaRes = fgseaRes[order(padj, -abs(NES)),]
  fgseaRes.save = fgseaRes
  
  
  for (i in 1: length(fgseaRes.save$leadingEdge)){
    collapesed = as.character(paste( unlist(fgseaRes.save$leadingEdge[i]), collapse=' - '))  
    fgseaRes.save$ENTREZ_GENES_IDS_IN[i] = collapesed
    
    
  }
  fgseaRes.save = fgseaRes.save[,-8]
  
  HELPER_SAVE_DATA_FRAME(as.data.frame(fgseaRes.save), paste0(outputDir,misegdbType,".csv"), rowNames = FALSE)
  
  return (list (result = fgseaRes, pathways = pathways))
  
}
ONTHO_GSEA_PlotPathWay = function (misegdbType, pathwayName, ranks,outputDir = NULL){
  if (!is.null(outputDir)){
    pdf (paste0(outputDir,misegdbType,"_", pathwayName,".pdf"))  
    plotEnrichment(misegdb[[pathwayName]], ranks)
    dev.off()
  }
  plotEnrichment(misegdb[[pathwayName]], ranks)
  
}
ONTHO_GSEA_Plot_Result = function (fgseaRes, padj = 0.05, title, outputDir = NULL, dbName ){
  #TODO Esto hay que agregarlo para visualizar de otra manera en fgssea
  fgseaResTidy <-fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj)
  
  
  
 
  graphics = ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=title ) + 
    theme_minimal()
  HELPER_SAVE_PDF(graphics, paste0(outputDir,dbName,".pdf"))
  
 
  
  fgseaRes.sig = fgseaResTidy[fgseaResTidy$padj<padj,]
  significants  = data.frame(fgseaRes.sig$NES)
  rownames(significants) = fgseaRes.sig$pathway
  
  pheatmap(significants ,cluster_rows = FALSE,cluster_cols = FALSE, filename = paste0(outputDir,dbName,"_heatmap.pdf"))
 
  
  
  
  #esta es otra manera
 # nUp = 20#length(fgseaRes$padj<padj & fgseaRes$ES >0)

  #nDown =20# length(fgseaRes$padj<padj & fgseaRes$ES <0)

  #topUp <- fgseaRes %>% 
   # filter(ES > 0) %>% 
    #top_n(nUp, wt=-padj)
  #topDown <- fgseaRes %>% 
   # filter(ES < 0) %>% 
    #top_n(nDown, wt=-padj)
  #topPathways <- bind_rows(topUp, topDown) %>% 
   # arrange(-ES)

  #if (!is.null(outputDir)){
   # pdf (paste0(outputDir,misegdbType,"_gseaTable.pdf"))  
    #plotGseaTable(pathways[topPathways$pathway], 
     #             ranks, 
      #            fgseaRes, 
       #           gseaParam = 0.5)
    
    #dev.off()
  #}
  
  #plotGseaTable(pathways[topPathways$pathway], 
   #             ranks, 
    #            fgseaRes, 
     #           gseaParam = 0.5)
  
  
}