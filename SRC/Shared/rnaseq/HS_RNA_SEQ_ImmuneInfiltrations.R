#TODO mergear la sample list para agregar el path a fpkm  y days to death
#cabmair not reported por NA |cambiar GROUP por sampleType | cambiar nombre de columnas " " por "_" y valores
#TODO Error in shapiro.test(z[, k]) : sample size must be between 3 and 5000
#TODO aca hay que poder elegir si FPKM o TPM y el path y que quiero corrrer dependiendo de donde vengo

#************************************************************************************************************************************************************************
# Load  FPKM o TPM
#************************************************************************************************************************************************************************
fpkmSampleList = HELPER_LoadCSV(samplesListFPKMPath)
fpkmSampleList = fpkmSampleList[fpkmSampleList$Sample_ID %in% samplesList$Sample_ID,]
fpkmSampleList = fpkmSampleList[order(match (fpkmSampleList$Sample_ID, samplesList$Sample_ID)),] #reordenamos por las dudas

htseqFPKM.rawCountMatrix = DESeq2_Merge_TCGA_ZippedExpressionFiles(sourceFolder = paste0(project.resourcePath,"rawData/HTSeq-FPKM"),
                                                               sampleFolder = fpkmSampleList$File_ID,
                                                               samplesFileName =fpkmSampleList$File_Name,
                                                               sampleID = fpkmSampleList$Sample_ID,
                                                               outputPath = NULL)

htseqTPM.rawCountMatrix = FPKMtoTPM(htseqFPKM.rawCountMatrix)

#*******************************************************************************************************
#merge annotation
#*******************************************************************************************************
geneAnnotationAttributes.immune = c('ensembl_gene_id','external_gene_name')
geneAnnotation.immune = ANN_GetAnnotationAll(specie = prject.rnaseq.ensemblSpecie, attributes = geneAnnotationAttributes.immune,collapseInRow = FALSE)

htseqFPKM.rawCountMatrix.annotation = HELPER_RowNamesAsFirstColumn(htseqFPKM.rawCountMatrix, "ENSENMBL_GENE_ID")
htseqFPKM.rawCountMatrix.annotation = ANN_MergeAnnotationToDataFrame(geneAnnotation.immune,htseqFPKM.rawCountMatrix.annotation,"ensembl_gene_id","ENSENMBL_GENE_ID", keepAll =  "Rigth")
htseqFPKM.rawCountMatrix.annotation = htseqFPKM.rawCountMatrix.annotation[,-1]

htseqTPM.rawCountMatrix.annotation = HELPER_RowNamesAsFirstColumn(htseqTPM.rawCountMatrix, "ENSENMBL_GENE_ID")
htseqTPM.rawCountMatrix.annotation = ANN_MergeAnnotationToDataFrame(geneAnnotation.immune,htseqTPM.rawCountMatrix.annotation,"ensembl_gene_id","ENSENMBL_GENE_ID", keepAll =  "Rigth")
htseqTPM.rawCountMatrix.annotation = htseqTPM.rawCountMatrix.annotation[,-1]
#*******************************************************************************************************

dir.create(file.path( paste0(project.resultsPath,"immune/Cibersort")), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path( paste0(project.resultsPath,"immune/xCell")), showWarnings = FALSE)
dir.create(file.path( paste0(project.resultsPath,"immune/quanTIseq")), showWarnings = FALSE)
dir.create(file.path( paste0(project.resultsPath,"immune/epic")), showWarnings = FALSE)

htseqFPKMPath = paste0(project.resultsPath,"immune/htseqFPKM.tsv")
htseqTPMPath = paste0(project.resultsPath,"immune/htseqTPM.tsv")
htseqFPKM.rawCountMatrix.annotation = htseqFPKM.rawCountMatrix.annotation[!is.na(htseqFPKM.rawCountMatrix.annotation$external_gene_name),]
htseqTPM.rawCountMatrix.annotation = htseqTPM.rawCountMatrix.annotation[!is.na(htseqTPM.rawCountMatrix.annotation$external_gene_name),]
HELPER_SaveTSV(htseqFPKM.rawCountMatrix.annotation,htseqFPKMPath)
HELPER_SaveTSV(htseqTPM.rawCountMatrix.annotation, htseqTPMPath)


#*******************************************************************************************************
#groups to compare for inmmune
#*******************************************************************************************************
project.groupsToCompare.inmune = project.groupsToCompare
for (i in 1:length(project.groupsToCompare.inmune)){
  project.groupsToCompare.inmune[[i]] = project.groupsToCompare.inmune[[i]][2:3] 
}
project.groupsToCompare.inmune
#*******************************************************************************************************
#Cibersort
#*******************************************************************************************************
source("./ExternalLibs/Cibersort/CIBERSORT.R")
immuneResults.Cibersort =  CIBERSORT("./ExternalLibs/Cibersort/LM22.txt",htseqFPKMPath, perm= 100, QN = FALSE, absolute=FALSE)




# immuneResults.Cibersort.filtered = immuneResults.Cibersort
# rownames(immuneResults.Cibersort.filtered) = gsub("\\.","-",rownames(immuneResults.Cibersort.filtered))
# 
# immuneResults.Cibersort.filtered = subset(immuneResults.Cibersort.filtered,rownames(immuneResults.Cibersort.filtered) %in% samplesList$Sample_ID)
# immuneResults.Cibersort.filtered = immuneResults.Cibersort.filtered[order(match (rownames(immuneResults.Cibersort.filtered),samplesList$Sample_ID)),]
# 
# immuneResults.Cibersort.result = IMMUNE_Humman_RNASEQ_STATS(immuneMatrix = immuneResults.Cibersort.filtered,
#                                                             pvalue = 0.05,
#                                                             outputDir =  paste0(project.resultsPath,"immune/Cibersort/") ,
#                                                             title =  "Cibersort TCGA - COAD",
#                                                             group = project.sample.group,
#                                                             decovolutionalMethod = "CIBERSORT"
# )
# 
# 
# setwd("/Volumes/Externo/Projects/PFI_ROADMAP/PRJ/")




immuneResults.Cibersort.result = IMMUNE_Humman_RNASEQ_STATS(immuneMatrix = immuneResults.Cibersort,
                                                            pvalue = 0.05,
                                                            outputDir = paste0(project.resultsPath,"./Cibersort/") ,
                                                            title =  "Cibersort TCGA - COAD",
                                                            group = project.sample.group,
                                                            decovolutionalMethod = "CIBERSORT",
                                                            groupToCompare = project.groupsToCompare.inmune,
                                                            orderFactor = project.groupsToCompare.inmune.orderFactor,
                                                            width = project.groupsToCompare.inmune.size.cibersort[1],
                                                            height = project.groupsToCompare.inmune.size.cibersort[2],
                                                            colCount = project.groupsToCompare.inmune.colCount.cibersort
                                                            
)





if (project.type == "TCGA"  & project.rnaseq.run.immune.survival){
  
  cibersort.survivalData = GET_IMMUNE_CIBERSORT_CLINICAL(immuneResults.Cibersort.result, samplesList)
  
  if (dim(cibersort.survivalData)[1]>1){
    GET_KAPLAN_MEIER_PLOT(survivalData = cibersort.survivalData,
                      outputDir = paste0(project.resultsPath, "immune/Cibersort/Survival/Kaplan-Meier/"))
    
    GET_COXPH_PLOT (survivalData = cibersort.survivalData,
                    outputDir =   paste0(project.resultsPath,"immune/Cibersort/Survival/coxph/") )
    
    # if (length(levels(project.sample.group)) >1){ #separamos tambien grupos a ver como da
    #   for (i in 1:length(levels(project.sample.group))){
    #     groupType = levels(project.sample.group)[i]
    #     groupNew = project.sample.group[samplesList$Sample_ID %in% rownames(cibersort.survivalData)]
    #     survivalData.filtered = cibersort.survivalData[groupNew == groupType, ]
    #     GET_KAPLAN_MEIER_PLOT(survivalData = survivalData.filtered,
    #                       outputDir = paste0(project.resultsPath, "immune/Cibersort/Survival/",groupType,"/"))
    #     
    #   }
    # }
  }
}



#***************************************************************
#xCell
#***************************************************************
devtools::install_github('dviraran/xCell')
library(xCell)
xCell.matrix = as.matrix(htseqFPKM.rawCountMatrix.annotation[,2:dim(htseqFPKM.rawCountMatrix.annotation)[2]])
rownames(xCell.matrix) = htseqFPKM.rawCountMatrix.annotation$external_gene_name
immuneResults.xCell = xCellAnalysis(xCell.matrix,rnaseq = TRUE)


immuneResults.xCell.result = IMMUNE_Humman_RNASEQ_STATS(immuneMatrix = t(immuneResults.xCell),
                                                            pvalue = 0.05,
                                                            outputDir =  paste0(project.resultsPath,"immune/xCell/") ,
                                                            title =  "xCell TCGA - COAD",
                                                            group = project.sample.group,
                                                            decovolutionalMethod = "XCELL",
                                                            groupToCompare = project.groupsToCompare.inmune,
                                                             orderFactor = project.groupsToCompare.inmune.orderFactor,
                                                            colCount =  project.groupsToCompare.inmune.colCount.xcell,
                                                             width = project.groupsToCompare.inmune.size.xcell[1],
                                                             height = project.groupsToCompare.inmune.size.xcell[2]
)

if (project.type == "TCGA"  & project.rnaseq.run.immune.survival){
  
  xCell.survivalData = GET_IMMUNE_CIBERSORT_CLINICAL(immuneResults.xCell.result, samplesList)
  
  if (dim(xCell.survivalData)[1]>1){
    GET_KAPLAN_MEIER_PLOT(survivalData = xCell.survivalData,
                      outputDir = paste0(project.resultsPath, "immune/xCell/Survival/Kaplan-Meier/"))

    
    GET_COXPH_PLOT (survivalData = xCell.survivalData,
                    outputDir =   paste0(project.resultsPath,"immune/xCell/Survival/coxph/") )
    
    
    # if (length(levels(project.sample.group)) >1){ #separamos tambien grupos a ver como da
    #   for (i in 1:length(levels(project.sample.group))){
    #     groupType = levels(project.sample.group)[i]
    #     groupNew = project.sample.group[samplesList$Sample_ID %in% rownames(xCell.survivalData)]
    #     survivalData.filtered = xCell.survivalData[groupNew == groupType, ]
    #     GET_KAPLAN_MEIER_PLOT(survivalData = survivalData.filtered,
    #                       outputDir = paste0(project.resultsPath, "immune/xCell/Survival/",groupType,"/"))
    #     
    #   }
    # }
  }
}
#***************************************************************
#quanTIseq  
#***************************************************************
htseqTPM.rawCountMatrix.annotation.unique = htseqTPM.rawCountMatrix.annotation
#filter duplicated symbol
htseqTPM.rawCountMatrix.annotation.unique = htseqTPM.rawCountMatrix.annotation.unique[!duplicated(htseqTPM.rawCountMatrix.annotation.unique$external_gene_name),]
HELPER_SaveTSV(htseqTPM.rawCountMatrix.annotation.unique, paste0(project.resultsPath,"immune/htseqTPMUnique.tsv"))

quanTIseq = paste0("sh ","./ExternalLibs/quanTIseq-master/quanTIseq_pipeline.sh  --inputfile=",project.resultsPath,"immune/htseqTPMUnique.tsv"," --pipelinestart=decon --tumor=TRUE --threads=7 --outputdir=",project.resultsPath,"immune/quanTIseq/")
quanTIseq
system(quanTIseq)

immuneResults.quantiseq = HELPER_LoadTSV(paste0(project.resultsPath,"immune/quanTIseq/quanTIseq_cell_fractions.txt"))
immuneResults.quantiseq = HELPER_FirstColumnAsRowName(immuneResults.quantiseq)


immuneResults.quantiseq.result = IMMUNE_Humman_RNASEQ_STATS(immuneMatrix = immuneResults.quantiseq,
                                                        pvalue = 0.05,
                                                        outputDir =  paste0(project.resultsPath,"immune/quantiseq/") ,
                                                        title =  "quantiseq TCGA - COAD",
                                                        group = project.sample.group,
                                                        decovolutionalMethod = "quantiseq",
                                                        groupToCompare = project.groupsToCompare.inmune,
                                                        orderFactor = project.groupsToCompare.inmune.orderFactor,
                                                        width = project.groupsToCompare.inmune.size.quantiseq[1],
                                                        height = project.groupsToCompare.inmune.size.quantiseq[2],
                                                        colCount = project.groupsToCompare.inmune.colCount.quantiseq
)

if (project.type == "TCGA"  & project.rnaseq.run.immune.survival){
  
  quanTIseq.survivalData = GET_IMMUNE_CIBERSORT_CLINICAL(immuneResults.quantiseq.result, samplesList)
  
  if (dim(quanTIseq.survivalData)[1]>1){
    GET_KAPLAN_MEIER_PLOT(survivalData = quanTIseq.survivalData,
                        outputDir = paste0(project.resultsPath, "immune/quanTIseq/Kaplan-Meier/"))
    
    
    GET_COXPH_PLOT (survivalData = quanTIseq.survivalData,
                    outputDir =   paste0(project.resultsPath,"immune/quanTIseq/Survival/coxph/") )
      # if (length(levels(project.sample.group)) >1){ #separamos tambien grupos a ver como da
      #   for (i in 1:length(levels(project.sample.group))){
      #     groupType = levels(project.sample.group)[i]
      #     groupNew = project.sample.group[samplesList$Sample_ID %in% rownames(quanTIseq.survivalData)]
      #     survivalData.filtered = quanTIseq.survivalData[groupNew == groupType, ]
      #     GET_KAPLAN_MEIER_PLOT(survivalData = survivalData.filtered,
      #                       outputDir = paste0(project.resultsPath, "immune/quanTIseq/Survival/",groupType,"/"))
      #     
      #   }
      # }
  }
}

#***************************************************************
#EPIC 
#***************************************************************
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
library(EPIC)
htseqTPM.rawCountMatrix.annotation.unique = htseqTPM.rawCountMatrix.annotation
#filter duplicated symbol
htseqTPM.rawCountMatrix.annotation.unique = htseqTPM.rawCountMatrix.annotation.unique[!duplicated(htseqTPM.rawCountMatrix.annotation.unique$external_gene_name),]
htseqTPM.rawCountMatrix.annotation.unique = HELPER_FirstColumnAsRowName(htseqTPM.rawCountMatrix.annotation.unique)

immuneResults.epic = EPIC(bulk = as.matrix(htseqTPM.rawCountMatrix.annotation.unique),TRef)


immuneResults.epic.result = IMMUNE_Humman_RNASEQ_STATS(immuneMatrix = immuneResults.epic$cellFractions,
                                                             pvalue = 0.05,
                                                             outputDir =  paste0(project.resultsPath,"immune/EPIC/") ,
                                                             title =  "EPIC TCGA - COAD",
                                                             group = project.sample.group,
                                                             decovolutionalMethod = "EPIC",
                                                       groupToCompare = project.groupsToCompare.inmune,
                                                       orderFactor = project.groupsToCompare.inmune.orderFactor,
                                                       width = project.groupsToCompare.inmune.size.epic[1],
                                                       height = project.groupsToCompare.inmune.size.epic[2],
                                                       colCount = project.groupsToCompare.inmune.colCount.epic
)

if (project.type == "TCGA"  & project.rnaseq.run.immune.survival){
  
  EPIC.survivalData = GET_IMMUNE_CIBERSORT_CLINICAL(immuneResults.epic.result, samplesList)
  
  if (dim(xCell.survivalData)[1]>1){
    GET_KAPLAN_MEIER_PLOT(survivalData = EPIC.survivalData,
                        outputDir = paste0(project.resultsPath, "immune/EPIC/Survival//Kaplan-Meier/"))
    
    
    GET_COXPH_PLOT (survivalData = EPIC.survivalData,
                    outputDir =   paste0(project.resultsPath,"immune/EPIC/Survival/coxph/") )
      # if (length(levels(project.sample.group)) >1){ #separamos tambien grupos a ver como da
      #   for (i in 1:length(levels(project.sample.group))){
      #     groupType = levels(project.sample.group)[i]
      #     groupNew = project.sample.group[samplesList$Sample_ID %in% rownames(EPIC.survivalData)]
      #     survivalData.filtered = EPIC.survivalData[groupNew == groupType, ]
      #     GET_KAPLAN_MEIER_PLOT(survivalData = survivalData.filtered,
      #                       outputDir = paste0(project.resultsPath, "immune/EPIC/Survival/",groupType,"/"))
      #     
      #   }
      # }
  }
}

#***************************************************************
#ESTIMATE  No funca bien
#***************************************************************
#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
#library(estimate)

#filterCommonGenes(input.f=output.counts.tpm.infiltrationMatrixFileName, output.f=output.infiltration.ESTIMATE.tempMatrix, id="GeneSymbol")
#estimateScore(output.infiltration.ESTIMATE.tempMatrix, output.infiltration.ESTIMATE.matrix, platform="illumina")
#immuneResults.estimate =loadTSV(output.infiltration.ESTIMATE.matrix)




