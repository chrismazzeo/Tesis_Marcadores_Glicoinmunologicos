#*******************************************************************************************************
#*******************************************************************************************************
#https://f1000research.com/articles/5-1384/v2
#*******************************************************************************************************
#*******************************************************************************************************
#TODO promerdiar las probesID para el mismo gen
#TODO Para GSE64658 no funca bien traer los grupos
#TODO que toma el db de la annotation automaticamente y no tener que especificar ni carglo en los settings
#TODO si es humano correr cibertsort
#TODO tsne.perplexity = 2, KMeans.kmax = 10,   KMeans.optimalNumberOfClusters = 4,  #TODO esto hay que automatizarlo
#TODO ver como traigo el platform para ver si es affy u otro e implementar los otros
#TODO  MAplot(esetObj,plot.method="smoothScatter", pairs=TRUE,groups=samplesList$GROUP)
#TODO falta agregar cibersort si es humano, hay otra cosa que no sea immucc?
#TODO ver como manejo de tener mas de un factor
#TODO ver bien si es por grupo si tengo que separar constrast por la correccion por multiple comparaciones
#TODO guardar el detalle del experimento
#TODO cuando devuelve celFiles tienen todo el path, hay que devolver solo el nombre de los archivos, y cuando leo affy pasarle el resourcePath
#GO y GSEA independiente de la especie
#TODO GSEA
#TODO DOSE
#*******************************************************************************************************
#*******************************************************************************************************


CreateInitialFolder = function (specie,projectName,groupsToCompare,resourcePath, resultPath ){
  dir.create(file.path( paste0(resourcePath,projectName)), showWarnings = FALSE)
  
  cellFolderPath = paste0( paste0(resourcePath,projectName),"/Cel/")
  dir.create(file.path(cellFolderPath), showWarnings = FALSE)
  
  dir.create(file.path( resultPath), showWarnings = FALSE)
  dir.create(file.path( paste0(resultPath,projectName)), showWarnings = FALSE)
  #results
  if (specie == "MOUSE"){
    dir.create(file.path( paste0(resultPath,projectName,"/ImmuCC/")), showWarnings = FALSE)  
  }
  dir.create(file.path( paste0(resultPath,projectName,"/Samples/")), showWarnings = FALSE)
  dir.create(file.path( paste0(resultPath,projectName,"/CountMatrix/")), showWarnings = FALSE)
  dir.create(file.path( paste0(resultPath,projectName,"/DescriptiveAnalysis/")), showWarnings = FALSE)
  dir.create(file.path( paste0(resultPath,projectName,"/DEA/")), showWarnings = FALSE)
  dir.create(file.path( paste0(resultPath,projectName,"/QA/")), showWarnings = FALSE)

  for (i in 1:length(groupsToCompare)){
    group = groupsToCompare[i]
    group = gsub("-", "_VS_", group) 
  
   
    dir.create(file.path( paste0(resultPath,projectName,"/DEA/",group)), showWarnings = FALSE) 
    dir.create(file.path( paste0(resultPath,projectName,"/DEA/",group,"/DEA")), showWarnings = FALSE) 
    dir.create(file.path( paste0(resultPath,projectName,"/DEA/",group,"/DEA/GlycoHeatmaps")), showWarnings = FALSE) 
    dir.create(file.path( paste0(resultPath,projectName,"/DEA/",group,"/DEA/GlycoHeatmaps/Log2FC_1")), showWarnings = FALSE) 
    
    
    dir.create(file.path( paste0(resultPath,projectName,"/DEA/",group,"/GO")), showWarnings = FALSE) 
    dir.create(file.path( paste0(resultPath,projectName,"/DEA/",group,"/GSEA")), showWarnings = FALSE) 
    dir.create(file.path( paste0(resultPath,projectName,"/DEA/",group,"/KEGG")), showWarnings = FALSE) 
  }
  return (paste0(resultPath,projectName,"/"))
}

Download_CEL = function (accessionNumber, source = "GEO", outputDir,projectName){
  #**********************************************************************
 
  celPath =  paste0(outputDir,projectName,"/Cel/")
  
  if (source == "GEO"){
      #1.Download RAW CEL FILE
      # Desde NCBI bajamos un Raw Cel File pasandole el Accession number
      # El Raw Cel file es un acrhivo comprimido que  contiene 1 Cel file por muestra 
      
      
      getGEOSuppFiles(accessionNumber,baseDir = celPath, makeDirectory = FALSE)
      
      #**********************************************************************
      filesToDelete = dir(celPath, full.names = TRUE)
      #1.1.unzip cel_raw file
      zippedCelPath= paste0(celPath,accessionNumber,"_raw.tar")
    
      untar(zippedCelPath, exdir= celPath)
      
      unlink(filesToDelete, recursive = FALSE)
      
      celFiles = dir(celPath, full.names = TRUE)
      return (celFiles)
       
  }
  else if(source == "ExpressArray") { #express array ENSEMBL
    #https://f1000research.com/articles/5-1384/v2
    anno_AE <- getAE(accessionNumber, path = cellPath, type = "raw")
    return (anno_AE)
  }
  else{  #local
    
  }

}
GEO_GET_SAMPLE_INFO = function(accessionNumber){
  gse <- getGEO(accessionNumber)
  gse = gse[[1]]
  group = gsub(" ", "_",gse@phenoData@data$source_name_ch1)
  group = gsub("-", "_",group)
  
  samplesList = data.frame(SAMPLE = gse@phenoData@data$description,ACCESSION_ID = gse@phenoData@data$geo_accession, FILES = "",DESCRIPTION = gse@phenoData@data$title , GROUP = group)
  return (samplesList)
}
AFFY_READ_CEL = function (files){
  #**********************************************************************
  #2.Read all CEL files 
  #arma un AfyBatch object que tienen los datos de las lecturas de las intensidades
  #tambien guardar otros datos que se pueden ver llamando a metodos ver ?AffyBatch
  #Podemos las probes, perfect match, miss match, el CDF
  
  #Se utiliza un chip determinado para cada especie y el resultado de expresion de genes se guarda en un Cel File
  #Probe (sondas) es un oligonucleotico de 25 bases que se aparea con un RNA target
  #la secuencias de las sondas se tomaron de GenBank, dbEST y RefSeq
  #Un chip contien hast 1.3 millones sondas de oligonueclotidos 
  #Probe cell, es cada area del chip donde se ubican las sondas, 
  
  #perfect match son las sondas que estan diseñaas para apareen perfectamente con el RNA target
  #PM es valor de la intesidad medida en las sondas donde hay perfect match
  #msmatch son las sondas que tienen una base distintas a la secuencia target y sirven para un apareanmento no especifico
  #MM es el valor de la intesidad medida en las donde donde hay mismatch
  #probe pair es una unida compuesta por una sonda perfect match y una sonda mismatch
  #probe  set contiene 11 probe pairs, que son los valores PMs y MMs relacionados con affyID
  #affyID es la identificacion para un probe set (que puede ser un gen o una fraccion del gen)
  # GeneChip Mouse Genome 430A 2.0 tienen 22600 probeset para analizar la expresion 14000 genes bien caracterizados
  # GeneChip Mouse Genome 430 2.0 tiene 45000 probe sets para analizar la expresion de 39000 transcriptos y sus variantes para mas de 34000 genes bien caracterisados
  # GeneChip Mouse Genome 430 2.0 puede correr hasta 96 muestras en paralelo
  
  
  #Cell File contiene las mediciones de intesidad e ubicacion para un micrroary donde fue hibridado

  #CDF   chip-definition file  contiene informacion del chip 
  #CDF file contiene la informacion relacionada con el probe pair set y su ubicacion The cdfName slot contains the necessary information for the package to find the locations of the probes for each probe set
  
  #varias sondas se  utilizan para represetar a los genes en la forma de probe set.
  #del Cel file tenemos para cada uno la ubicacion fisica 
  #el CEL file tambien contiene el nombre del CDF file nque se usa para mapear la ubicacion de los probe set
  #el CDF guarda los probe set relacionados a cada ubicacion en chip
  
  #The ProbeSet class holds the information of all the probes related to an affyID. The components are pm and mm.
  #The method probeset extracts probe sets from AffyBatch objects. It takes as arguments an AffyBatch object and a vector of affyIDs and returns a list of objects of class ProbeSet
  
  
  
  #Leemos todos los cel files de la carpeta y nos devuelve un objeto del tipo AffyBatch 
  #AffyBatch class:
  #expresion:  una matriz de n muestras x 1004004 mediciones
  #cantidad de probes en el chip,
  #cantidad de muestras
  #cantidad de genes que detecta el chip
  #cdf cdfName(affydata)
  #Annotation  son los affyID
  

  affydata <- affy::ReadAffy(filenames = files) #depende la marca, aca usamos  microarray de Affymetrix

  return (affydata)
}

AFFY_GET_EXPRESSION_MATRIX = function (affyData, backgroundCorretionType,outputDir){
  #*******************************************************************************************************
  #Expression matrix
  #*******************************************************************************************************
  #Crea un objeto eset de la clase ExpressionSet que hereda de eSet  
  #Que hace
  #1. lee loss probe level data  del objeto affyData y los convierte a expression values, generando una matriz (gene AffyID-sample)
  #2. background correction  RMA o GCRMA
  #3 normalization
  #4. probe specific background correction ej remueve los MM
  #5 sumariza los probe set values en un valor de expresion
  
  #background correction usadno RMA 
  #1. Probe specific correction of the PM probes using a model based on observed intensity being the sum of signal and noise
  #2. Normalization of corrected PM probes using quantile normalization (Bolstad et al.,2003)
  #. Calculation of Expression measure using median polish.
  # En el chip  se incluyen un set genes para facilitar la normalizacion y la esclamiento,  que sirven para despues hacer comapraciones entre las lecturas de expresion 
  # entre distintos chips. Estos genes tienen una expresion constante en diversos tejidos (son housekeeping)
  
  
  if (backgroundCorretionType == "GCRMA"){
    esetObj <- gcrma(affyData)
  }else{
    esetObj <- affy::rma(affyData)  
  }

  #Guardamos la expression matrix orginal
  esetObjSave = as.data.frame(exprs(esetObj))
  esetObjSave = HELPER_RowNamesAsFirstColumn(esetObjSave, "PROBE_ID")
  HELPER_SAVE_DATA_FRAME(esetObjSave,paste0(outputDir,"CountMatrix/CountMatrix.csv"))

  return (esetObj)
}

Affy_QualityControl = function (affyData,expressionObj,outputDir){
  #4. Quality Control
 # image(affyData)
  
  #Verificamos los valores de expresion antes y despues de la normalizacion para que las muestras sean comparables
  pdf(paste0(outputDir,"QA/Normalization.pdf"))
  cols <- brewer.pal(8, "Set1")
  par(mfrow = c(2, 2))
  boxplot(affyData, col=cols, main = "Unnormalised intensity values", las = 2)
  # plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles.gcrma
  boxplot(expressionObj, col=cols, main = "Normalised intensity values",las = 2)
  
  # Plot a density vs log intensity histogram for the unnormalised data
  hist(affyData, col=cols,main = "Unnormalised intensity values")
  # Plot a density vs log intensity histogram for the normalised data
  hist(expressionObj, col=cols ,main = "Normalised intensity values")
  dev.off()

  #RNA degradation QA
  par(mfrow = c(1,1))
  #se ordena los probes dentro de un probeset por ubicacion relativa al extremo 5' del RNa target
  #Como la degracion del RNA tipicamente comienza desde el extremo 5' de la molecula, esperariamos que las intesidades de las sondas sean siempre mas bajas en el extremo 5' que en el exteremo 3'
  deg <- AffyRNAdeg(affyData)
  summaryAffyRNAdeg(deg)
  pdf(paste0(outputDir,"QA/RNADegradation.pdf"))
  plotAffyRNAdeg(deg)
  dev.off()
 
  
  
  pdf(paste0(outputDir,"QA/ProbeIntensity.pdf"))
  #expressed genes by filtering based on intensity --> agregar AffyQuality --> esto se puede usar para remover baja intensidad pero hay que poner un cutoff a mano
  medians <- rowMedians(Biobase::exprs(expressionObj))
  hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                   main = "Histogram of the median intensities",
                   border = "antiquewhite4",
                   xlab = "Median intensities")
  dev.off()
}

AFFY_FILTER_LOW_INTENSITY = function(expressionObj,samplesGroup, minValue = 0,outputDir){
  #removemos los probeset con baja intensidad para que no haga probelmas
  if (minValue >0){
 
    no_of_samples <- table(samplesGroup)
    samples_cutoff <- min(no_of_samples)

    
    idx_man_threshold <- apply(Biobase::exprs(esetObj), 1,
                               function(x){
                                 sum(x > minValue) >= samples_cutoff})
    table(idx_man_threshold)
    expressionObj <- subset(expressionObj, idx_man_threshold)
    pdf(paste0(outputDir,"QA/ProbeIntensityFiltered.pdf"))
    #expressed genes by filtering based on intensity --> agregar AffyQuality --> esto se puede usar para remover baja intensidad pero hay que poner un cutoff a mano
    medians <- rowMedians(Biobase::exprs(expressionObj))
    hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                     main = "Histogram of the median intensities",
                     border = "antiquewhite4",
                     xlab = "Median intensities")
    dev.off()  
    return (expressionObj)
    
  }
  
}
AFFY_GET_IMMUCC_MATRIX = function(esetObj,annotationDB, sampleName,sampleGroup, outputDir){ 
  #*******
  # El resultado lo tengo que subir a la web http://218.4.234.74:3200/immune/#
  #******
  #The microarray platform further involves Affymetric mouse 430 2.0, Illumina MouseWG-6 v2.0 expression beadchip and Agilent Whole Mouse Genome Microarray 4x44K v2.  


  library(affy)                      # Version: 1.56.0
  library(frma)                      # Version: 1.30.1
  library(mouse4302mmentrezgcdf)     # Version: 19.0.0
  library(mouse4302frmavecs)         # Version: 1.5.0
  library(preprocessCore)            # Version: 1.40.0
  library(mouse4302mmentrezg.db)
  # Note: the library version listed here is the package versions used in my work. However, you can choose the correspondent version based on the R installed in your computer and it will have no impact on the result.
  
  # Read all cel files under path with a custom cdf "mouse4302mmentrezgcdf"
  affydata <- ReadAffy(filenames= celFiles, cdfname = "mouse4302mmentrezgcdf")
  
  # Preprocessing with frma
  eset <- frma(affydata)
  
  # Output the expression value of samples profiled on array
  ematrix <- exprs(eset)
  
  annotations <- AnnotationDbi::select(
    x       = mouse4302mmentrezg.db,
    keys    = featureNames(eset),
    columns = c("ENTREZID","SYMBOL"), #columns = c("PROBEID","ENTREZID", "SYMBOL"),
    keytype = "PROBEID"
  )
  
  matrix.entrez <- data.frame(ENTREZID = annotations$ENTREZID,ematrix, stringsAsFactors = FALSE)
  matrix.entrez = matrix.entrez[!is.na(matrix.entrez$ENTREZID),]
  table(duplicated(matrix.entrez$ENTREZID))
  HELPER_SAVE_DATA_FRAME(matrix.entrez,paste0(outputDir,"/immucc/immuccMatrixSource.tsv"))

}

MICROARRAY_GET_ANNOTATION = function(annotationDbName, annotationDB){
  #annotations
  #Hay 2 problemas 
  #1. Hay probesID duplicados que mapean con distintos genes --> los saco
  #2 distintos probesID mapean con el mismo gen --> por el sesgo de la posicion del probeset o porque hibridan con splicing alternativo, o porque estan diseñadas para hibridar con distintos exones del transcripto --> que hago?
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1784106/
  #solucion--> ordeno por B y me quedo con el que tenga mayor valor de B --> es un parametro que indica si que tanto esta diferencialmente expresado
  
  #columns(project.microarray.annotationDB)
  
  #//////// annotations
 
  annotationDbName = paste0(annotationDbName,".db")
  #BiocManager::install(annotationDbName, version = "3.8", update = FALSE,ask = FALSE)
  library(package = annotationDbName,character.only = TRUE)

  annotation = select(project.microarray.annotationDB, keys = keys(project.microarray.annotationDB), columns = c("SYMBOL"))
  table(duplicated(annotation$PROBEID))
  
  #removemos los probesID que matchean para 2 o mas genes distintos
  duplicatedProbes = annotation[duplicated(annotation$PROBEID),]
  
  finalAnnotation = annotation[!(annotation$PROBEID %in% duplicatedProbes$PROBEID),]
  table(duplicated(finalAnnotation$PROBEID))
  finalAnnotation = finalAnnotation[!is.na(finalAnnotation$SYMBOL),]
  dim(finalAnnotation)
  
  #guardamos los probesId que tenian multimapping solo por las dudas
  duplicatedProbes <- group_by(annotation, PROBEID)
  duplicatedProbes = dplyr::summarize(duplicatedProbes, no_of_matches = n_distinct(SYMBOL))
  duplicatedProbes = duplicatedProbes[order(duplicatedProbes$no_of_matches,decreasing = TRUE),]
  #HELPER_SAVE_DATA_FRAME(annotationDuplicatedSumarized, paste0(outputDir,"CountMatrix/DuplicatedProbesID.csv"))
  return (finalAnnotation)
  
}

MICROARRAY_PLOT_RESULTS = function(result,padj,groupToCompare,  outputDir ){
 # DESEQ2_volcano(result,padj,filePath = paste0(outputDir,"DEA/Volcano_plot.pdf"))
  DESEQ2_log2FCdistribution(result,padj, outputDir = paste0(outputDir,"DEA/"), fileName = "LOG2Distribution.pdf")
  
  
  title = project.groupsToCompare[i]
  title = gsub("-", " VS ", title) 
  pdf(paste0(outputDir,"DEA/pvalueHistogram.pdf"))
  hist(result$pvalue, col = brewer.pal(3, name = "Set2")[1],
       main = title, xlab = "p-values")
  dev.off()
  
 
}



project.resultsPath = CreateInitialFolder(specie = project.specie,
                                          projectName = project.name,
                                          groupsToCompare = project.groupsToCompare,
                                          resourcePath = project.resourcePath ,
                                          resultPath = project.resultsPath
)

geneFilteringList =  HELPER_LoadCSV(geneListFilterPath)
HELPER_SAVE_DATA_FRAME(geneFilteringList, paste0(project.resultsPath,"GeneFilteringList.csv")) 




celFiles = Download_CEL(accessionNumber = project.microarray.downloadAccessionNumber, #download CEL files
                        source = project.microarray.cel.source,
                        outputDir = project.resourcePath,
                        projectName = project.name)

if (!is.na(samplesListPath)){
  samplesList = HELPER_LoadCSV(samplesListPath)

}else{
  samplesList = GEO_GET_SAMPLE_INFO(project.microarray.downloadAccessionNumber)
  samplesList$FILES = celFiles
}
samplesList$GROUP = factor(samplesList$GROUP)
HELPER_SAVE_DATA_FRAME(samplesList, paste0(project.resultsPath,"/samples/SampleList.csv"), rowNames = FALSE)

affyData = AFFY_READ_CEL(files = celFiles) # read cel files


esetObj = AFFY_GET_EXPRESSION_MATRIX(affyData = affyData,
                                     backgroundCorretionType =  project.microarray.backgroundCorrectionMode,
                                     outputDir = project.resultsPath) # background correction, normalize, summarization

Affy_QualityControl(affyData = affyData,
                    expressionObj = esetObj,
                    outputDir = project.resultsPath)  #Completar y guardar graficos

DA_MICROARRAY_ALL(rawCountMatrix = esetObj,
                  samplesNames = samplesList$SAMPLE,
                  samplesGroups = samplesList$GROUP,
                  tsne.perplexity = 2,  #TODO esto hay que automatizarlo
                  KMeans.kmax = dim(esetObj)[2]/2,  
                  KMeans.optimalNumberOfClusters = length(levels(factor(samplesList$GROUP))), 
                  KMeans.maxIteration = settings.kmeans.numberOfinteractions,
                  KMeans.numbersOfRandomCenters = settings.kmeans.numbersOfRandomCenters,
                  groups.colorArray = settings.color,
                  outputDir = project.resultsPath)

esetObj = AFFY_FILTER_LOW_INTENSITY (expressionObj = esetObj,
                                     samplesGroup = samplesList$GROUP,
                                     minValue = project.microarray.settings.lowIntensityFilter ,
                                     outputDir =  project.resultsPath)

annotation = MICROARRAY_GET_ANNOTATION (annotationDbName =  annotation(affyData), 
                                        annotationDB = project.microarray.annotationDB)
#******** IMMUCC ***********************
if (!is.na(project.immucc.svrPath)){
  
  immuCCResult.svr = IMMUNE_ImmuCC_RNASEQ_STATS(filePath = project.immucc.svrPath,
                                                outputDir = paste0(project.resultsPath,"immuCC/") ,
                                                title =  "ImmuCC SVR",
                                                samplesName = samplesList$SAMPLE,
                                                group = samplesList$GROUP,
                                                decovolutionalMethod = "SVR"
  )
}
if (!is.na(project.immucc.llsrPath)){
  
  immuCCResult.svr = IMMUNE_ImmuCC_RNASEQ_STATS(filePath = project.immucc.llsrPath,
                                                outputDir = paste0(project.resultsPath,"immuCC/") ,
                                                title =  "ImmuCC LLSR",
                                                samplesName = samplesList$SAMPLE,
                                                group = samplesList$GROUP,
                                                decovolutionalMethod = "LLSR"
  )
}


#**************COMPARACIÖN ***********
#asigno samples a grupos, una pelotudez esto y armo una matriz
design<- model.matrix(~0+GROUP, data = samplesList) #ver si hay que agregar mas grupos
colnames(design) =gsub("GROUP", "", colnames(design)) #removemos el Group que se no se porque se creo


#Calculo del promedio de cada gen de cada grupo
fit = lmFit(esetObj, design) 

#grupos que voy a comparar
contrast.matrix = makeContrasts(contrasts = project.groupsToCompare,
                                levels = samplesList$GROUP)

fit.contrast = contrasts.fit(fit,contrast.matrix)

#Test estadisticos
fit.eb = eBayes(fit.contrast) #me devuelve todas las comparaciones de los grupos



#correcion por multiple comparaciones por cada comparacion, guardamos y downstream analysis
for (i in 1:dim(fit.eb)[2]){
  print (project.groupsToCompare[i])
  outputDirGroup = paste0(project.resultsPath,"DEA/",project.groupsToCompare[i],"/")
  outputDirGroup = gsub("-", "_VS_", outputDirGroup) 
  
  
  tempResults = topTable(fit.eb,coef = i,adjust.method="BH",number = Inf, lfc = project.settings.DEA.H0.Log2FC)#
  tempResults = HELPER_RowNamesAsFirstColumn(tempResults,"PROBEID")
  
  colnames(tempResults)[2] = "log2FoldChange"
  colnames(tempResults)[5] = "pvalue"
  colnames(tempResults)[6] = "padj"
  
  #saco todos los probeID que no estan anotados o estaban duplicados en las anotaciones
  tempResults.unique = tempResults[tempResults$PROBEID %in% annotation$PROBEID,]
  tempResults.annotation = ANN_MergeAnnotationToDataFrame(annotation,tempResults.unique,"PROBEID","PROBEID", keepAll = "Left")
  head(tempResults.annotation)
  
  
  
  #filtro probesIDs que mapean al mismo gen y me quedo con el más significativo
  tempResults.annotation = tempResults.annotation[order(tempResults.annotation$B, decreasing = TRUE),]
  tempResults.annotation.final = tempResults.annotation[!duplicated(tempResults.annotation$SYMBOL),]
  
  
  colnames(tempResults.annotation.final)[2] = "external_gene_name"
  head(tempResults.annotation.final)
  
  #cambiar el nombre log2FoldChange  external_gene_name
  
  
  HELPER_FILTER_SAVE_PFI(tempResults.annotation.final,
                         geneList = geneFilteringList, 
                         padj = project.settings.DEA.padj, 
                         log2FC = project.settings.DEA.cutoff.log2FC, 
                         outputDir = outputDirGroup,
                         project.type = project.type)
  
  
  MICROARRAY_PLOT_RESULTS (result = tempResults.annotation.final,
                           padj = project.settings.DEA.padj,
                           groupToCompare = project.groupsToCompare[i],
                           outputDir = outputDirGroup)
  
  #*********************
  #GO
  
  DE_GENES_ID = tempResults$PROBEID[tempResults$padj<project.settings.DEA.padj] #todos las probeID  diferencialmentes expresados
  
  #esto es solo para graficar 
  #aca buascamos el background para ver que genes se expresan maso menos igual (mas trucho es hacer como background todos los que no se expresan diferencialemente)
  back_genes_idx <- genefilter::genefinder(exprs(esetObj),
                                           as.character(DE_GENES_ID),
                                           method = "manhattan", scale = "none")
  back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
  back_genes <- featureNames(esetObj)[back_genes_idx]
  back_genes <- setdiff(back_genes, DE_GENES_ID)#termino teniendo la lsita de genes que maso menos se expresan igual menos los que se expresasn diferencialmente, asi obtengo el background
  
  pdf(paste0(outputDirGroup,"GO/backgroundVsDE.pdf"), width = 10, height = 20)
  multidensity(list(
    all = tempResults[,"AveExpr"] ,
    fore = tempResults[DE_GENES_ID , "AveExpr"],
    back = tempResults[rownames(tempResults) %in% back_genes, "AveExpr"]),
    col = c("#e46981", "#ae7ee2", "#a7ad4a"),
    xlab = "mean expression",
    main = "DE genes for CD-background-matching")
  dev.off()
  
  #TOP GO
  geneList = ONTHO_GO_GENE_LIST(tempResults,padj = project.settings.DEA.padj, log2FC_UP = 1, log2FC_DOWN = 1)
  
  top_GO_data <- new("topGOdata", ontology = "BP", allGenes = geneList,
                     nodeSize = 10, annot = annFUN.db, affyLib = paste0(annotation(affyData),".db"))
  
  #The algorithm starts processing the nodes / GO categories from the highest (bottommost) level and then iteratively moves to nodes from a lower level. If a node is scored as significant, all of its genes are marked as removed in all ancestor nodes. This way, the “elim” algorithm aims at finding the most specific node for every gene.
  result_top_GO_elim =     runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
  
  result_top_GO_classic =  runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")
  
  res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                         Fisher.classic = result_top_GO_classic,
                         orderBy = "Fisher.elim" , topNodes = 100)
  
  #summary
  genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                             chip = paste0(annotation(affyData),".db"), geneCutOff = 1000)
  
  
  #Los genes que estan en una term y song significativos tienen  pvalue = 2
  res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
    str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"]," - "),
          collapse = "")
  })
  res_top_GO$Term=gsub(",", " - ", res_top_GO$Term)
  HELPER_SAVE_DATA_FRAME(res_top_GO, paste0(outputDirGroup,"GO/EnrichedTerms.csv"))
  pdf(paste0(outputDirGroup,"GO/EnrichedTerms.pdf"), width = 10, height = 20)
  graphicsOver =  res_top_GO %>% 
    top_n(n = 100, wt=-`Rank in Fisher.classic`) %>% 
    mutate(hitsPerc=Significant*100/Annotated) %>% 
    ggplot(aes(x=hitsPerc, 
               y=Term, 
               colour=`Rank in Fisher.classic`, 
               size=Significant)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="GO term", colour="Rank", size="Count")

  print(graphicsOver)
  dev.off()
  
#  showSigOfNodes(top_GO_data, score(result_top_GO_elim), firstSigNodes = 100,
  #               useInfo = 'def')
  
  #DISEASE ONTHOLOGY   DOSE
  #mapeo a ENTREZID
 # entrez_ids <- mapIds(project.microarray.annotationDB,
  #                     keys = names(geneList),
   #                    keytype = "PROBEID",
    #                   column = "ENTREZID")
  
  #universe = entrez_ids
  #universe = as.vector(universe[!is.na(names(universe)) & !is.na(universe)])
  
  #deGenes = entrez_ids[DE_GENES_ID]
  #deGenes = as.vector(deGenes[!is.na(names(deGenes)) & !is.na(deGenes)])
  
  
  #reactome_enrich <- enrichPathway(gene = deGenes,
   #                                universe = universe,
    #                               organism = "mouse",
     #                              pvalueCutoff = 0.05
  #)
  #reactome_enrich@result$Description <- paste0(str_sub(
   # reactome_enrich@result$Description, 1, 20),
  #  "...")
  
  #head(as.data.frame(reactome_enrich))[1:6]
  #barplot(reactome_enrich)
  #emapplot(reactome_enrich, showCategory = 10)
}
unlink(paste0(project.resourcePath,project.name,"/Cel"), recursive = TRUE)
print ("******* FINISHED ******** ")
