

#******** MERGE expression files **********

#********************************
# From GDC, download samples files 
#
#Input:
#   sourceFolder = folder containing htseq files
#   Samples = the data.frame associated to the files (fileName, case id, etc)
#
#Output:
#      data.frame with rows = genesId and columns samples gene expression
#
#Note  columns names = sampleID, if some sampleID is not unique, it will add to the sampleID the postFix "_Rep"
#********************************
DESeq2_Merge_TCGA_ZippedExpressionFiles = function (sourceFolder,sampleFolder,samplesFileName,sampleID, outputPath = NULL){
  
  #filesInFolder = dir(sourceFolder,full.names = TRUE, recursive = TRUE)
  #filesToCopy = filesInFolder[grep(".gz",filesInFolder)]
  
  
  HTCountsTempFolder = paste0(sourceFolder,"/tempFolder")
  dir.create(path = HTCountsTempFolder)
  
  filesToCopy = paste0(sourceFolder,"/",sampleFolder,"/",samplesFileName)
  file.copy(from = filesToCopy, to = HTCountsTempFolder)
  gunzipCommand = paste0("gunzip ", HTCountsTempFolder,"/* --force")
  system(gunzipCommand)
  
  #saco el .gz final
  fileName = unlist(strsplit(paste0(HTCountsTempFolder,"/",samplesFileName),"\\.gz",0))
  
  #armo el count matrix
  outputDataFrame = read_tsv(fileName[1], col_names = FALSE)
  outputDataFrame = as.data.frame(outputDataFrame)
  colnames(outputDataFrame) =c("Ensembl_Gene_ID",as.character(sampleID[1]))
  
  for (i in 2:length(fileName)){
    dat = read_tsv(fileName[i], col_names = FALSE)
    dat = as.data.frame(dat)
    dat = dat[-1]
    colnames(dat) = as.character(sampleID[i])
    outputDataFrame = cbind(outputDataFrame,dat)
  }
  outputDataFrame = outputDataFrame[order(outputDataFrame$Ensembl_Gene_ID,decreasing = FALSE),]
  outputDataFrame = DESeq2_RemoveHTSeqExtraInfo(outputDataFrame)
  
  outputDataFrame$Ensembl_Gene_ID = DESeq2_RemoveGeneVersion(outputDataFrame$Ensembl_Gene_ID)
  rownames(outputDataFrame) = outputDataFrame$Ensembl_Gene_ID
  outputDataFrame = outputDataFrame[,-1]
  
  HELPER_SAVE_DATA_FRAME(outputDataFrame,outputPath)
  
  unlink(HTCountsTempFolder, recursive = TRUE)
  
  return (outputDataFrame)
}

#para outputs de htseq individuales sin compactar
DESeq2_MergeHTSeqExpressionFiles = function(sourceFolder,fileNames, outputPath = NULL){
  filePath = paste0(sourceFolder,"/",fileNames)
  #armo el count matrix
  outputDataFrame = read_tsv(filePath[1], col_names = FALSE)
  outputDataFrame = as.data.frame(outputDataFrame)

  colnames(outputDataFrame) =c("Ensembl_Gene_ID",as.character(fileNames[1]))

  for (i in 2:length(filePath)){
    dat = read_tsv(filePath[i], col_names = FALSE)
    dat = as.data.frame(dat)
    dat = dat[-1]
    colnames(dat) = as.character(fileNames[i])
    outputDataFrame = cbind(outputDataFrame,dat)
  }
  outputDataFrame = outputDataFrame[order(outputDataFrame$Ensembl_Gene_ID,decreasing = FALSE),]
  outputDataFrame = DESeq2_RemoveHTSeqExtraInfo(outputDataFrame)
  
  outputDataFrame$Ensembl_Gene_ID = DESeq2_RemoveGeneVersion(outputDataFrame$Ensembl_Gene_ID)
  rownames(outputDataFrame) = outputDataFrame$Ensembl_Gene_ID
  outputDataFrame = outputDataFrame[,-1]
  #  outputDataFrame = outputDataFrame[,-1]
  
  HELPER_SAVE_DATA_FRAME(outputDataFrame,outputPath)
  return (outputDataFrame)
}
#*************
#VST 
#vemos la varianza de cada gen para cada muestra - el promedio para ese gen
#*************

#******* DESEQ2 *********

# rawExpressionDataFrame = rowsnames(geneID), colnames(SampleIDs), rows x Columns = raw expression
# conditionDataFrame  rowsnames(sampleID), columns = factor
# ncol (rawExpressionDataFrame) == nrow (conditionDataFrame)
# formula = ~ condition  (choose factors or factors)
# coresNumber  for parallel processing >1
#prefilteringCuttOff a value >0

DESeq2_tx_import = function(filesPath, transcriptAnnotation, type = "kallisto", geneLevel, conditionDataFrame,formula, coresNumber = 1, prefilteringCuttOff = 0){
  #hay 2 tipos de archivos por cada sample, abundance.h5 y abundance.tsv
  #abundance.tsv tiene 4 filas target_id	length	eff_length	est_counts	tpm
  #abundance.h5 tiene 4 dataframe y el bootstrap
  #abundance  --> es el tpm del tsv
  #counts --> es el est_counts del tsv
  #length  --> el el length del tsv
  #countsFromAbundance
  
  
  #type recibe que metodo se uso para contar los reads
  #txOut = false sumariza los transcriptos a  nivel gen  | si usamos O TRUE luego podemos sumarizar usando  txi.sum <- summarizeToGene(txi.tx, tx2gene)
  
  
  txi<- tximport(filesPath,tx2gene = transcriptAnnotation, type = type, txOut = !geneLevel)
  dds = DESeqDataSetFromTximport(txi = txi, 
                                 colData = conditionDataFrame,
                                 design = formula)
  
  #prefiltering
  if (prefilteringCuttOff >0){
    nrow(dds)
    dds <- dds[ rowSums(counts(dds)) > prefilteringCuttOff, ]
    nrow(dds)
  }
  
  if (coresNumber >1){
    register(MulticoreParam(coresNumber))
  }
  dds <- DESeq(dds, parallel = (coresNumber >1))
  return (dds)
}
DESeq2_Matrix_Import = function (rawExpressionDataFrame, conditionDataFrame, formula,  coresNumber = 1, prefilteringCuttOff = 0){
  
  #creamos la matriz con todo
  dds = DESeqDataSetFromMatrix (countData = rawExpressionDataFrame,
                                colData = conditionDataFrame,
                                design = formula)
  

  #prefiltering
  if (prefilteringCuttOff >0){
    nrow(dds)
    dds <- dds[ rowSums(counts(dds)) > prefilteringCuttOff, ]
    nrow(dds)
  }

  #***************************************************************
  #DF analysis
  #***************************************************************
  #utilizo paralalelismo
  if (coresNumber >1){
    register(MulticoreParam(coresNumber))
  }
  #hago el analisis de expresion diferencial
  #This very simple function call does all the hard work. Briefly, this function performs three things:
  #estimation of size factor(Compute a scaling factor for each sample to account for differences in read depth and complexity between samples
  #the estimation of dispersion values for each gene
  #Test for differences in expression among groups 
  
  #independent filtering, lo que hace es remover los genes con low count que son muy dificiles de ver si  hay expresion diferencial, pero si 
  #los dejamos va a afectar en la correcion por multiple testing, entonces al removerlos antes de hacer la correccion por FDR, vamos a tener
  #comor esultado mas genes significativos
  #deseq2 lo hace automaticamente
  
  dds <- DESeq(dds, parallel = (coresNumber >1))
  resultsNames(dds)

  return (dds)
}

DESeq2_VST = function(dds){
  #variance stabilization
  #es una normalizacion de deseq para hacer analsis descriptivo
  vsd <- vst(dds, blind = FALSE)
}

# coresNumber  for parallel processing >1
# adjusthedMethod =  see ?p.adjust
#lfcThreshold = if you want to change the null hipothesis for a LFC
#contrast = c ("var factor", "factoA", "factorB")
#filePath path to sabe results
#filter fun  #Independent hypothesis weighting , supuedamente es mejor que BH, si le pasamos uno, no hace BH

#NOTA
# greaterABS vs greater, en greaterABS es 2 tails, esto es que prueba el log2fc tanto + como -, mientras que en greater, solo es el positivo, por eso hay que usar greaterABS
# greaterABS = greater + less

DESeq2_ExtractResults = function(dds, filterFun = NULL, coresNumber = 1, pAdjustMethod = "BH", adjustedPValue = 0.05, contrast = NULL, lfcThreshold = 0, format ="DataFrame"){
  if (is.null(contrast)){
    print ("You need to specify the contast")
    return (NULL)
  }
  #***************************************************************
  #Comparo de a 2 grupos
  #***************************************************************
  #compare 2 groups
  #seteamos:
  #método de correcion por comparacion multiple
  # alpha
  # que grupos comparamos
  if (!is.null(filterFun)){
    res <- results(dds, parallel = (coresNumber >1), pAdjustMethod = pAdjustMethod, contrast=contrast, lfcThreshold = lfcThreshold, format = format,alpha = adjustedPValue, filterFun = filterFun)
  }
  else{
    res <- results(dds, parallel = (coresNumber >1), pAdjustMethod = pAdjustMethod, contrast=contrast, lfcThreshold = lfcThreshold, format = format,alpha = adjustedPValue)
  }
  #Nota
  # si tengo mas grupos para comaprar res.2 =  y armo el otro contrast
  
  #El resultado es una matriz con:
  #1. promedio
  #2. log2FoldChange  = log2(treated / untreated)  |||| 2^foldchange
  #3. pValue 
  #4. adjusted pValue (es lo que usamos para rechazar H0)

  
 
  
  return (res)
}

DESeq2_RemoveGeneVersion = function(data){
  data = gsub("\\..*","",data)
  return (data)
}
DESeq2_RemoveHTSeqExtraInfo = function(data){
  
  #remove last 5 rows 
  if ("__no_feature" %in% data$Ensembl_Gene_ID){
    data = data[-c(1:5),]
  }
  return (data)
}


DESeq2_Result_ALL_Stats = function (DESeqResult,DESeqResult.shrunken = NULL, DESeqResult.annotation = NULL, padj = 0.05, log2FCCutOff = 1,outputDir){

  DESEQ2_summary(DESeqResult,filePath = paste0(outputDir,"DEA/summary.txt"))
  
  DESEQ2_plotMA(DESeqResult,padj, filePath = paste0(outputDir,"DEA/MA_Plot.pdf"))
  
  if (!is.null(DESeqResult.shrunken)){
    #Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes
    #It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary 
    #filtering thresholds.
   
    DESEQ2_plotMA(DESeqResult.shrunken,padj, filePath = paste0(outputDir,"DEA/MA_Plot_Shrunken.pdf"))
  }
  DESEQ2_volcano(DESeqResult,padj,filePath = paste0(outputDir,"DEA/Volcano_plot.pdf"))
  
  DESEQ2_log2FCdistribution(DESeqResult,padj, outputDir = paste0(outputDir,"DEA/"), fileName = "LOG2Distribution.pdf")
  
  if (!is.null(DESeqResult.annotation)){
     DESEQ2_biotypeDistribution(DESeqResult.annotation, outputDir = paste0(outputDir,"DEA/"), fileName = "BioTypeDistribution.pdf")
  }
  
 # DESEQ2_DE_Chromosome(DESeqResult.annotation, filePath = paste0(outputDir,"DEA/chromosome.csv"))
}

DESEQ2_DE_Chromosome= function(DESeqResult.annotation,filePath){

  total = as.data.frame(table(DESeqResult.annotation$chromosome_name))
  colnames(total) = c("CHROMOSOME", "Genes")
  total = total[-(20:42),]
  total
  
  totalDe = as.data.frame(table(DESeqResult.annotation$chromosome_name[DESeqResult.annotation$padj<0.5]))
  if (dim(totalDe)[1]>0){
    colnames(totalDe) = c("CHROMOSOME", "DE")
    totalDe = totalDe[-(20:29),]
    totalDe
    final = merge(total,totalDe, by = "CHROMOSOME")
    final = data.frame(final, GeneRatio = round(final$DE/final$Genes,1))
   
    HELPER_SAVE_DATA_FRAME(final,filePath)
  }
  print(final)
}
DESEQ2_biotypeDistribution = function(DESeqResult.annotation, outputDir, fileName){
  biotype = as.data.frame(table(DESeqResult.annotation$gene_biotype))
  if(dim(biotype)[1] >0){
    graphics = plot_ly(biotype, type = "bar", x=~Var1, y = ~Freq)%>%
      layout(title = "Biotype Distribution",
             xaxis = list(title = "",
                          zeroline = TRUE),
             yaxis = list(title = "Freq")
    )
    HELPER_SAVE_ORCA (graphics, outputDir = outputDir, fileName = fileName)
  }

}
DESEQ2_log2FCdistribution = function(DESeqResult, padj,outputDir,fileName){
  #log2fc distribution
  all = DESeqResult$log2FoldChange[!is.na(DESeqResult$padj) & DESeqResult$padj<=padj]
  all = as.data.frame(table(round(all,0)))
  print(all)
  if(dim(all)[1] >0){
    graphics = plot_ly(all, type = "bar", x=~Var1, y = ~Freq)%>%
      layout(title = "Log2FC distribution of significant genes",
             xaxis = list(title = "Log2FC",
                          zeroline = TRUE),
             yaxis = list(title = "Freq")
      )
  
    HELPER_SAVE_ORCA (graphics, outputDir = outputDir, fileName = fileName)
  }

}
DESEQ2_summary = function(DESeqResult,filePath = NULL){
  #***************************************************************
  #summary results
  #***************************************************************
  #describmos la informacion sobre los test que usamos
  #baseMean average of the normalized count values, divided by the size factors, taken over all samples in the DESeqDataSet
  #log2FoldChange is the effect size estimate. It tells us how much the gene’s expression seems to have changed due to treatmen
  #lfcSE es la estimacion del estandar error asociado al log2 fold change
  #The purpose of a test for differential expression is to test whether the data provides sufficient evidence to conclude that this
  #value is really different from zero. DESeq2 performs for each gene a hypothesis test to see whether evidence is sufficient to decide 
  #against the null hypothesis that there is zero effect of the treatment on the gene and that the observed difference between treatment
  #and control was merely caused by experimental variability
  #pvalue result os the test, que tan posible que el fold change para un gen sea lo que indica el fold change
  #padj es el pvalue corregido por multiples comparaciones  y es el que realmente importa
  
  #nota si el pvalue = NA es porque el gene tenia todo 0 o porque habia muchos outliers
  if (!is.null(filePath)){
    sink(filePath)
    print ("***************************************************************")
    print ("DE Analysis Results")
    print ("***************************************************************")
    
    mcols(DESeqResult)$description[2]
    print ("***************************************************************")
    print ("Summary")
    
    summary(DESeqResult)
    sink(type = "message")
    sink()
  }
  
  print (mcols(DESeqResult)$description[2])
  print(summary(DESeqResult))
}

DESEQ2_plotMA = function(DESeqResult,padj, filePath){
  #Plot MA
  #In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. 
  #Points will be colored red if the adjusted p value is less than alpha1. Points which fall out of the window are plotted as open triangles pointing either up or down.
  pdf(filePath)
  DESeq2::plotMA(DESeqResult,alpha = padj, main = "MA plot")
  abline(1,0, col="dodgerblue", lwd=2)
  abline(-1,0, col="dodgerblue", lwd=2)
  abline(2,0, col="dodgerblue", lwd=2,lty=2)
  abline(-2,0, col="dodgerblue", lwd=2,lty=2)
  axis(side=2, at=c(-10:10))
  dev.off()
  
  DESeq2::plotMA(DESeqResult,alpha = padj, main = "MA plot")
  abline(1,0, col="dodgerblue", lwd=2)
  abline(-1,0, col="dodgerblue", lwd=2)
  abline(2,0, col="dodgerblue", lwd=2,lty=2)
  abline(-2,0, col="dodgerblue", lwd=2,lty=2)
  
  axis(side=2, at=c(-10:10))
}
DESEQ2_volcano = function(DESeqResult,cutoff = 0.05 ,title = "Volcano Plot", minLog2FC = 1,colors = c("gray","yellow","blue","red"), xMaxValues = c(10,10), yMaxValues = c(0,10),showLegend = TRUE,geneLabelMinLog2FC = 2, geneLabelMinAlpha = 0.05,showlines = FALSE, showGeneNames = TRUE,showGeneNameList = NA,labelOffset = 0.5,drawLabelBackground = false,filePath){
  
  #Nota sobre pValues = NA puede ser debido a 
  #1. If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
  #2. If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. 
  #3. If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. 
  
  #Volcano plot

 
  pdf(filePath)
  filePath
  

#  abline(h=-log10(cutoff), col="black",lty=4, lwd=2.0)
  DESeqResult$log10adjPvalues = -log10(DESeqResult$padj)
  DESeqResult$log10adjPvalues[-log10(DESeqResult$padj)>yMaxValues[2]] = yMaxValues[2]
 
  
  #Gray no significant genes
  with(DESeqResult, plot(log2FoldChange, log10adjPvalues, cex=0.6,pch=20, main=title, xlab="Log2 fold change", ylab="-log10(adjusted p-value)",xlim=xMaxValues,ylim=yMaxValues,col = colors[1]))
  if (showlines){
    abline(v=c(-1,1), col="black", lty=2)
    abline(h=-log10(cutoff), col="black", lty=2)
  }
  
  #axis(side=1, at=c(-log2FC:log2FC))
  #axis(side=2, at=c(-100))
 
  # Yellow up/down genes with Log2FC <|1|
  with(subset(DESeqResult, padj<cutoff &(log2FoldChange > -minLog2FC & log2FoldChange < minLog2FC)), points(log2FoldChange, log10adjPvalues, cex=0.6,pch=20, col= colors[2]))
  
  #Blue down genes  with Log2FC <|1|
  with(subset(DESeqResult, padj<cutoff & log2FoldChange <= -minLog2FC), points(log2FoldChange, log10adjPvalues, cex=0.6,pch=20, col= colors[3]))
  
  #Red up genes  with Log2FC >|1|
  with(subset(DESeqResult, padj<cutoff & log2FoldChange >= minLog2FC), points(log2FoldChange, log10adjPvalues, cex=0.6,pch=20, col= colors[4]))

  #labels
  if (showGeneNames){
    #up
    
    if (all(is.na(showGeneNameList))){
      gn.selectedUp <- DESeqResult$log2FoldChange > geneLabelMinLog2FC & DESeqResult$padj < geneLabelMinAlpha   
      gn.selectedDown <- DESeqResult$log2FoldChange < -geneLabelMinLog2FC & DESeqResult$padj < geneLabelMinAlpha 
    }
    else{
      DESeqResult = DESeqResult[DESeqResult$external_gene_name %in% showGeneNameList,]
      gn.selectedUp <- DESeqResult$log2FoldChange > geneLabelMinLog2FC & DESeqResult$padj < geneLabelMinAlpha 
      gn.selectedDown <- DESeqResult$log2FoldChange < -geneLabelMinLog2FC & DESeqResult$padj < geneLabelMinAlpha 
    }
  
    if (drawLabelBackground){
      #up
      for (i in 1:length(gn.selectedUp)){
        if (gn.selectedUp[i]){
          legend(DESeqResult$log2FoldChange[i] ,
                 -log10(DESeqResult$padj)[i] + labelOffset, DESeqResult$external_gene_name[i], box.col = "black", 
                 bg = "white", adj = 0.4,cex = 0.7)
        }
      }
      #down
      for (i in 1:length(gn.selectedDown)){
        if (gn.selectedDown[i]){
          legend(DESeqResult$log2FoldChange[i] ,
                 -log10(DESeqResult$padj)[i]+ labelOffset, DESeqResult$external_gene_name[i], box.col = "white", bg = "white", adj = 0.4,cex = 0.7)
        }
      }
    }
    else{
      
      #up
      text(DESeqResult$log2FoldChange[gn.selectedUp],
         -log10(DESeqResult$padj)[gn.selectedUp] + labelOffset,
         lab=DESeqResult$external_gene_name[gn.selectedUp ],
         cex=0.8)
    
      #down
    
      text(DESeqResult$log2FoldChange[gn.selectedDown] ,
        -log10(DESeqResult$padj)[gn.selectedDown] + labelOffset,
        lab=DESeqResult$external_gene_name[gn.selectedDown ],
       cex=0.8)
    }
  
    
    
  }

  #legend
  if (showLegend){
    legend("topright", legend=c( "padj<0.05 & |Log2FC| < 1", "Down-regulated","Up-regulated"),
         fill=c( "yellow", "blue","red"), horiz=FALSE, cex=0.8)  
  }
  dev.off()

}
sDESEQ2_volcano2 = function(DESeqResult,cutoff = 0.05 ,title = "Volcano Plot", minLog2FC = 1,color = "gray", xMaxValues = c(10,10), yMaxValues = c(0,10),showLegend = TRUE,geneLabelMinLog2FC = 2, geneLabelMinAlpha = 0.05,showGeneNames = TRUE,filePath){
  
  pdf(filePath)
  DESeqResult$log10adjPvalues = -log10(DESeqResult$padj)
  DESeqResult$log10adjPvalues[-log10(DESeqResult$padj)>yMaxValues[2]] = yMaxValues[2]
  
   cols <- densCols(DESeqResult$log2FoldChange, DESeqResult$log10adjPvalues)
   plot(DESeqResult$log2FoldChange, DESeqResult$log10adjPvalues, col=cols, panel.first=grid(),
        main=title, xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
        pch=20, cex=0.6,xlim=xMaxValues,ylim=yMaxValues)
   abline(v=0)
   abline(v=c(-1,1), col="brown")
   abline(h=-log10(cutoff), col="brown")
   
   #labels
   if (showGeneNames){
     gn.selected <- abs(DESeqResult$log2FoldChange) > geneLabelMinLog2FC & DESeqResult$padj < geneLabelMinAlpha 
     text(DESeqResult$log2FoldChange[gn.selected],
          -log10(DESeqResult$padj)[gn.selected],
          lab=DESeqResult$external_gene_name[gn.selected ], cex=0.6)
   }
   #legend
   if (showLegend){
     legend("topright", legend=c( "padj<0.05 & |Log2FC| < 1", "Down-regulated","Up-regulated"),
            fill=c( "yellow", "blue","red"), horiz=FALSE, cex=0.8)  
   }
   dev.off()
 }
# #***************************************************************
# #Session information
# #***************************************************************
# sink(output.packageFileName)
# print ("***************************************************************")
# print ("Package utilizados")
# 
# devtools::session_info()
# sink(type = "message")
# sink()

#***************************************************************
#Removing hidden batch effects con SVA
#***************************************************************
#es <- results(dds, contrast=c("dex","trt","untrt"))
# if(removeBathEffect){
#   
#   folder.results.batchEffect = paste0(folder.results.alpha,"batchEffect/")
#   dir.create(file.path(folder.results.batchEffect), showWarnings = FALSE)
#   
#   
#   
#   library("sva")
#   dat  <- counts(dds, normalized = TRUE)
#   
#   #filter low counts
#   idx  <- rowMeans(dat) > 1
#   dat  <- dat[idx, ]
#   #full model with gene expression
#   fullModell  <- model.matrix(~ condition, data =colData(dds))
#   #null model aca es para ver las adjustment variables ej pheno data
#   nullModel <- model.matrix(~   1, colData(dds))
#   #no pasamos n.sv para que lo calcule automaticamente
#   svseq <- svaseq(dat, fullModell, nullModel)
#   #sv es una mtrix donde las collumnas corresponde a la estimacion de variables subrrogadas
#   #pprob.gam es la probabilitadad posterior de cada gene este asociadio a 1 o mas variables latente
#   #pprob.b es la problabildiad posterior de cada gene esta asociado a las variables de interes
#   #n.sv es el numero de surrogate variables que nos devuelve sva, y es lo que usamos para agregar a nuetra modelo
#   
#   print (paste0("number of subrrogate variables: ",svseq$n.sv))
#   #svseq$sv
#   par(mfrow = c(svseq$n.sv, 1), mar = c(3,5,3,1))
#   
#   pdf(paste0(folder.results.batchEffect,"batchEffectSurrogateVariables.pdf"))
#   for (i in 1:svseq$n.sv) {
#     stripchart(svseq$sv[, i] ~ dds$condition, vertical = TRUE, main = paste0("SV", i))
#     abline(h = 0)
#   }
#   dev.off()
#   #esto es solo para analysis
#   #ahora vemos expression differencial comparamo el full model vs el null model sin agregar ninguna variable subrogada
#   pValues = f.pvalue(dat,fullModell,nullModel)
#   qValues = p.adjust(pValues,method="BH")
#   table(qValues<alphaValue)
#   
#   #ahora hacemos lo mismo pero agregando las variables subrogadas
#   modSv = cbind(fullModell,svseq$sv)
#   mod0Sv = cbind(nullModel,svseq$sv)
#   pValuesSv = f.pvalue(dat,modSv,mod0Sv)
#   qValuesSv = p.adjust(pValuesSv,method="BH")
#   table(qValuesSv<alphaValue)
#   
#   #ahora corro DESeq2 con las variables subrrogadas
#   ddssva <- dds
#   ddssva$SV1 <- svseq$sv[,1]
#   ddssva$SV2 <- svseq$sv[,2]
#   design(ddssva) <- ~ SV1 + SV2 + condition
#   
#   #direct adjunstmen
#   #si conozco el batch effect ej para gtext-tcga el source seria un batch effect, aca puedo aplicar directgamente combat
#   batch = pheno$batch
#   modcombat = model.matrix(~1, data=colData(dds))
#   combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE)
#   
#   pValuesComBat = f.pvalue(combat_edata,mod,mod0)
#   qValuesComBat = p.adjust(pValuesComBat,method="BH")
#   
#   #Removing known batch effects with a lin- ear model
#   modBatch = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
#   mod0Batch = model.matrix(~as.factor(batch),data=pheno)
#   pValuesBatch = f.pvalue(edata,modBatch,mod0Batch)
#   qValuesBatch = p.adjust(pValuesBatch,method="BH")
#   
#   #***************************************************************
#   #Removing hidden batch effects con RUVSeq
#   #***************************************************************
#   library("RUVSeq")
# }
# 
# 

#DEXSEQ
DEXSEQ_CALCULATE = function (files, sampleList, flattenedfile, design,fitExpToVar, coresNumber = 1){
 
  
  dexseq.res =  DEXSeqDataSetFromHTSeq(
    countfiles=files,
    sampleData = as.data.frame(samplesList),
    design = design,
    flattenedfile=flattenedfile )
  if (coresNumber >1){
    BPPARAM = MulticoreParam(workers=coresNumber)
  }
  DEXSeq(dexseq.res,BPPARAM=BPPARAM, fitExpToVar = "GROUP")
  
  dexseq.res = DEXSeq(dexseq.res, BPPARAM=BPPARAM, fitExpToVar =fitExpToVar)
}
DEXSEQ_PLOT_ALL_HTML = function(dexseq.res,padj = 0.05, fitExpToVar, outputDir, fileName,coresNumber = 1,specie){
  if (coresNumber >1){
    BPPARAM = MulticoreParam(workers=coresNumber)
  }
  
  mart = useEnsembl(biomart="ensembl", dataset=specie)
  
  DEXSeqHTML(object = dexseq.res,fitExpToVar = fitExpToVar, path = outputDir, file = fileName, FDR = padj,BPPARAM = BPPARAM, mart = mart, attributes = attributes)
  
}
DEXSEQ_PLOT_GENE = function (dexseq.res, geneID, fitExpToVar,padj = 0.05, outputDir = NULL){

  if (!is.null(outputDir)){

      pdf(paste0(outputDir,geneID,".pdf"),width = 12,height = 12)
      par(mar=c(1,1,1,1))
      plotDEXSeq( dexseq.res, geneID = geneID,fitExpToVar = fitExpToVar,  cex.axis=1.2, cex=1.3, lwd=1,FDR = padj, norCounts = TRUE, expression = TRUE,splicing = TRUE, displayTranscripts = TRUE , names = TRUE, legend = TRUE)
      dev.off()
      
  }
  else{
     print(plotDEXSeq( dexseq.res, geneID = geneID,fitExpToVar = fitExpToVar,  cex.axis=1.2, cex=1.3, lwd=1,FDR = padj, norCounts = TRUE, expression = TRUE,splicing = TRUE, displayTranscripts = TRUE , names = TRUE, legend = TRUE))
  }
}
DEXSEQ_PLOT_ALL_PDF = function(dexseq.res,listID, fitExpToVar, padj,outputDir){
  for (i in 1:length(listID)){
    try(
      DEXSEQ_PLOT_GENE(dexseq.res, geneID = listID[i], fitExpToVar= fitExpToVar,padj = padj, outputDir = outputDir)
    )
  }
  
  
}

