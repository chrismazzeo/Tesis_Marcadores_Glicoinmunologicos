#clear console
cat("\014") 
#clear enviroment
#rm(list=ls())
list.of.packages <- c("RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(RColorBrewer)
#******************************************************************************************************************************
#Project Settings /// hay que seleccionar desde una bd
#******************************************************************************************************************************
Sys.setenv('MAPBOX_TOKEN' = 'pk.eyJ1IjoiY2hyaXNtYXp6ZW8iLCJhIjoiY2puZ2NsM3dxMDFycTNwbzUwa3ZrdXp1byJ9.epVfBOfjT3GSkZG5OqsXQA')
project.source = "GDC"
project.program = "TCGA"
project.subprogram = "COAD"

project.folder = paste0(project.source,"_",project.program,"_",project.subprogram)
folder.input.resources = paste0("./resources/",project.source,"/",project.program,"/",project.subprogram,"/")
#******************************************************************************************************************************
#Input files Settings
#******************************************************************************************************************************
input.geneMouseFilteringFileName = "./Resources/ListaDeGenesMM.csv"
input.geneHSFilteringFileName = "./Resources/ListaDeGenesHS.csv"
input.geneHS_MMFilteringFileName = "./Resources/ListaDeGenesHS-MM.csv"
input.clinical =  paste0(folder.input.resources,"/clinical/","clinical.tsv")
input.clinical.exposure =  paste0(folder.input.resources,"/clinical/","exposure.tsv")
input.clinical.msistatus=  paste0(folder.input.resources,"/clinical/","All_CDEs.txt")
input.clinical.therapies =  paste0(folder.input.resources,"/clinical/","clinical_drug_public_coad.txt.xls")

input.reads.htSeq.countsFile = paste0(folder.input.resources,"Reads/","HTSeq-Counts-Samples.tsv")
input.reads.htSeq.FPKMFile = paste0(folder.input.resources,"Reads/","HTSeq-FPKM-Samples.tsv")

#******************************************************************************************************************************
#Input folders
#******************************************************************************************************************************
folder.input.reads.htSeq.countsFile = paste0(folder.input.resources,"Reads/","HTSeq-Counts/")
folder.input.reads.htSeq.FPKMFile = paste0(folder.input.resources,"Reads/","HTSeq-FPKM/")

#******************************************************************************************************************************
#Folder Output Settings
#******************************************************************************************************************************
folder.output.results = ("./Results/")
dir.create(file.path(folder.output.results), showWarnings = FALSE)

folder.output.results = paste0(folder.output.results,project.folder,"/")
dir.create(file.path(folder.output.results), showWarnings = FALSE)

#****************************
#Samples
#****************************
folder.output.Samples = paste0(folder.output.results,"Samples/")
dir.create(file.path(folder.output.Samples), showWarnings = FALSE)

#****************************
#descriptive analysis
#****************************
folder.output.da = paste0(folder.output.results,"DescriptiveAnalysis/")
dir.create(file.path(folder.output.da), showWarnings = FALSE)

#****************************
#Count Matrix
#****************************
folder.output.countMatrix = paste0(folder.output.results,"CountMatrix/")
dir.create(file.path(folder.output.countMatrix), showWarnings = FALSE)

#****************************
#Differential expression
#****************************
folder.output.de = paste0(folder.output.results,"DE/")
dir.create(file.path(folder.output.de), showWarnings = FALSE)

#****************************
#Infiltration
#****************************
folder.output.infiltration = paste0(folder.output.results,"Infiltrations/")
dir.create(file.path(folder.output.infiltration), showWarnings = FALSE)

folder.output.infiltration.cibersort = paste0(folder.output.infiltration,"Cibersort/")
dir.create(file.path(folder.output.infiltration.cibersort), showWarnings = FALSE)

folder.output.infiltration.xcell = paste0(folder.output.infiltration,"XCELL/")
dir.create(file.path(folder.output.infiltration.xcell), showWarnings = FALSE)

folder.output.infiltration.quanTIseq = paste0(folder.output.infiltration,"quanTIseq/")
dir.create(file.path(folder.output.infiltration.quanTIseq), showWarnings = FALSE)

folder.output.infiltration.timer = paste0(folder.output.infiltration,"TIMER/")
dir.create(file.path(folder.output.infiltration.timer), showWarnings = FALSE)

# folder.output.infiltration.epic = paste0(folder.output.infiltration,"EPIC/")
# dir.create(file.path(folder.output.infiltration.epic), showWarnings = FALSE)

# folder.output.infiltration.estimate = paste0(folder.output.infiltration,"ESTIMATE/")
# dir.create(file.path(folder.output.infiltration.estimate), showWarnings = FALSE)


folder.output.infiltration.bioespecimen = paste0(folder.output.infiltration,"Bioespecimen/")
dir.create(file.path(folder.output.infiltration.bioespecimen), showWarnings = FALSE)

#****************************
#Heatmaps
#****************************
folder.output.heatmaps = paste0(folder.output.results,"Heatmaps/")
dir.create(file.path(folder.output.heatmaps), showWarnings = FALSE)

folder.output.heatmapsByFamilyGroup = paste0(folder.output.heatmaps,project.folder,"FamilyGroup/") 
dir.create(file.path(folder.output.heatmapsByFamilyGroup), showWarnings = FALSE)

#****************************
#Annotations
#****************************
folder.output.annotations = paste0(folder.output.results,"Annotations/")
dir.create(file.path(folder.output.annotations), showWarnings = FALSE)


#******************************************************************************************************************************
#Output files settings
#******************************************************************************************************************************
#annotations
#****************************
output.annotations.hs = paste0(folder.output.annotations,"HomoSapiensGeneAnnotation.tsv")
output.annotations.mouse = paste0(folder.output.annotations,"MouseGeneAnnotationOrthologs.tsv")
#****************************
#Samples Filtered
#****************************
output.htseqCountsSamples = paste0(folder.output.Samples,"HTSeqRawCountsSamples.tsv")
output.htseqFPKMSamples  = paste0(folder.output.Samples,"HTseqFPKMSamples.tsv")
output.clinicalSamples  = paste0(folder.output.Samples,"clinicalSamples.tsv")
#*************************
# Expression Matrix
#*************************
#countMatrix
output.countMatrix.FileName = paste0(folder.output.countMatrix,"RawCountMatrix.csv")
output.countMatrix.FilteredBySignificantFileName =  paste0(folder.output.countMatrix,project.folder,"_CountMatrix(FilteredBy_SignificantGenes)")

#FPKM Matrix
output.countMatrix.FPKM.FileName = paste0(folder.output.countMatrix,"FPKMCountMatrix.csv")
#TPM Matrix
output.countMatrix.TPM.FileName = paste0(folder.output.countMatrix,"TPMCountMatrix.csv")
#*************************
#DE Results
#*************************
output.DE.FileName =  paste0(folder.output.de,project.folder,"_DE_Results.csv")
output.DE.FilteredBySignificantFileName =paste0(folder.output.de,project.folder,"_DE_Results(FilteredBy_SignificantGenes).csv")
output.DE.Glyco.HS.FileName =paste0(folder.output.de,project.folder,"_DE_Results(Glyco_HS).csv")
output.DE.Glyco.HS_MM.FileName =paste0(folder.output.de,project.folder,"_DE_Results(Glyco_HS_MM).csv")
#*************************
#infiltration matrix
#*************************
output.counts.raw.infiltrationMatrixFileName = paste0(folder.output.countMatrix,"RawCounts_InfiltrationMatrix.csv")
output.counts.fpkm.infiltrationMatrixFileName = paste0(folder.output.countMatrix,"FPKM_InfiltrationMatrix.csv")
output.counts.tpm.infiltrationMatrixFileName = paste0(folder.output.countMatrix,"TPM_InfiltrationMatrix.csv")


output.infiltration.cibersort = paste0(folder.output.infiltration.cibersort,"CibersortResult.csv")
output.infiltration.xcell = paste0(folder.output.infiltration.xcell,"xCellResult.csv")

output.infiltration.quanTiSeq.temp = paste0(folder.output.infiltration.quanTIseq,"quanTIseq_cell_fractions.txt")
output.infiltration.quanTiSeq = paste0(folder.output.infiltration.quanTIseq,"quanTIseqResult.csv")

output.infiltration.TIMER = paste0(folder.output.infiltration.timer,"TIMER.csv")

#output.infiltration.epic = paste0(folder.output.infiltration.epic,"EPIC.csv")

# output.infiltration.ESTIMATE.tempMatrix = paste0(folder.output.infiltration.estimate,"Temp_ESTIMATE_Matrix.tsv")
# output.infiltration.ESTIMATE.matrix = paste0(folder.output.infiltration.estimate,"ESTIMATE Matrix.tsv")
# output.infiltration.ESTIMATE.plot = paste0(folder.output.infiltration.estimate,"ESTIMATE Plot_",project.folder,".pdf")

#*************************
#file description
output.filesDescriptionFileName = paste0(folder.output.results,"README_FIRST.txt")
#*************************

#*************************
#package description
output.packageFileName = paste0(folder.output.results,"packageUtilizados.txt")
#*************************
#******************************************************************************************************************************
#Color Settings
#******************************************************************************************************************************
settings.color = c(  
  
  "chartreuse1","deeppink3","cyan3","yellow", "darkblue","darkorchid1","darkgoldenrod1",
  "chartreuse4","hotpink","cadetblue1","darkorchid4","coral3", "brown4","azure4","darkgoldenrod4",
  "forestgreen","deeppink4","cadetblue3","chocolate1","coral1","brown3","cornsilk4",
  "chartreuse3", "deeppink1","cyan1","darkviolet", "orangered", "chocolate3","azure3","blue",
  
  "darkolivegreen1","hotpink4","aquamarine3","brown1","darkgoldenrod3",
  "darkolivegreen3","cyan4","aquamarine4",  
  "darkolivegreen4", "cadetblue4","chocolate4","cornsilk1","black",
  "chartreuse1","deeppink3","cyan3","yellow", "darkblue","darkorchid1","darkgoldenrod1",
  "chartreuse4","hotpink","cadetblue1","darkorchid4","coral3", "brown4","azure4","darkgoldenrod4",
  "forestgreen","deeppink4","cadetblue3","chocolate1","coral1","brown3","cornsilk4",
  "chartreuse3", "deeppink1","cyan1","darkviolet", "orangered", "chocolate3","azure3","blue",
  
  "darkolivegreen1","hotpink4","aquamarine3","brown1","darkgoldenrod3",
  "darkolivegreen3","cyan4","aquamarine4",  
  "darkolivegreen4", "cadetblue4","chocolate4","cornsilk1","black",
  "chartreuse1","deeppink3","cyan3","yellow", "darkblue","darkorchid1","darkgoldenrod1",
  "chartreuse4","hotpink","cadetblue1","darkorchid4","coral3", "brown4","azure4","darkgoldenrod4",
  "forestgreen","deeppink4","cadetblue3","chocolate1","coral1","brown3","cornsilk4",
  "chartreuse3", "deeppink1","cyan1","darkviolet", "orangered", "chocolate3","azure3","blue",
  
  "darkolivegreen1","hotpink4","aquamarine3","brown1","darkgoldenrod3",
  "darkolivegreen3","cyan4","aquamarine4",  
  "darkolivegreen4", "cadetblue4","chocolate4","cornsilk1","black"
  
)
#******************************************************************************************************************************
#PCA Settings
#******************************************************************************************************************************
output.pca.2dFileName = paste0(folder.output.da,project.folder,"_PCA_2D.pdf")
output.pca.2dWithClustersFileName = paste0(folder.output.da,project.folder,"_PCA_2DWithClusters.pdf")
output.pca.3dFileName = paste0(folder.output.da,project.folder,"_PCA_3D.png")
output.pca.3dWithClustersFileName = paste0(folder.output.da,project.folder,"_PCA_3DWithClusters.png")
output.pca.PCPorcentage = paste0(folder.output.da,project.folder,"_PCPorcentage.pdf")

#******************************************************************************************************************************
#kmeans Settings
settings.kmeans.maxCluster = 3
settings.kmeans.maxIteration = 100


output.kmeans.optimalClustersFileName = paste0(folder.output.da,project.folder,"_Kmeans_optimalClusters.pdf")
output.kmeans.resultFileName =  paste0(folder.output.da,project.folder,"_kmeans.pdf")
output.kmeans.custerFileName = paste0(folder.output.da,project.folder,"_Clusters.csv")



#******************************************************************************************************************************
#TSNE
#******************************************************************************************************************************
output.tsne.2d.resultFileName =  paste0(folder.output.da,project.folder,"_tsne2d.pdf")
output.tsne.2d.optimalClusterFileName =  paste0(folder.output.da,project.folder,"_tsne2dOptimalCluster.pdf")
output.tsne.3d.resultFileName = paste0(folder.output.da,project.folder,"_tsne3d.pdf")
output.tsne.3d.optimalClusterFileName =  paste0(folder.output.da,project.folder,"_tsne3dOptimalCluster.pdf")

#******************************************************************************************************************************
#Differential expression analysis settings
#******************************************************************************************************************************
settings.de.prefilteringCuttOff = 0
settings.de.alphaValue = 0.05
settings.de.lfcThreshold = 1
settings.de.summaryFileName = paste0(folder.output.de,"desummary.txt")

#******************************************************************************************************************************
#Correlation
#******************************************************************************************************************************
settings.correlation.cutoff = 0.6
settings.correlation.normalityResultFileName = paste0(folder.output.da,"normalityResult.txt")
#******************************************************************************************************************************
#heatmap Settings
#******************************************************************************************************************************
settings.heatmap.maxColumn = 130
settings.heatmap.naColor = "#black"
settings.heatmap.de.color=c("blue","black","red")
settings.heatmap.correlation.color=c("green","white","brown")
settings.heatmp.distance.colorRamp = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
settings.heatmap.cellWidth = 6
settings.heatmap.cellHeight = 3

output.heatmap.sampleDistanceFileName = paste0(folder.output.da,"SAMPLE_DISTANCE_HEATMAP.pdf")

output.heatmap.vst.GeneFilteredSifnificantExpressionFileName = paste0(folder.output.heatmaps,project.folder,"_VST_ExpressionAcrossSamples(FilterdByGeneList_Significant)_HEATMAP.pdf")

output.heatmap.countMatrix.CorrelationGeneFilteredCorrelationFileName  = paste0(folder.output.heatmaps,project.folder,"_countMatrix_Correlation(FilterdByGeneList_Significant)_HEATMAP.pdf")
output.heatmap.countMatrix.CorrelationGeneFilteredCorrelationHighCorrelationFileName  = paste0(folder.output.heatmaps,project.folder,"_countMatrix_Correlation(FilterdByGeneList_Significant_HighCorrelation)_HEATMAP.pdf")
output.heatmap.countMatrix.CorrelationGeneFilteredCorrelationHighCorrelationSummaryFileName  = paste0(folder.output.heatmaps,project.folder,"_countMatrix_CorrelationSummary(FilterdByGeneList_Significant_HighCorrelation)_HEATMAP.pdf")

output.heatmap.deResultsFileName = paste0(folder.output.heatmaps,project.folder,"_DE_RESULTS_HEATMAP.pdf")
output.heatmap.deSignificantsFileName  = paste0(folder.output.heatmaps,project.folder,"_DE_SignificantGenes_HEATMAP.pdf")
output.heatmap.deSignificantsFilteredFileName  = paste0(folder.output.heatmaps,project.folder,"(FilteredBy_GeneList)_HEATMAP.pdf")
output.heatmap.deSignificantsFilteredLFCFileName  = paste0(folder.output.heatmaps,project.folder,"(FilteredBy_GeneList__FoldChangeCutOff)_HEATMAP.pdf")
output.heatmap.deSignificantsFilteredLFCFileName  = paste0(folder.output.heatmaps,project.folder,"AsList(FilteredBy_GeneList__FoldChangeCutOff)_HEATMAP.pdf")

