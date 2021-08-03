source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "TCGA"
project.specie ="HUMAN"
project.name = "COAD_SampleType"   
project.resourcePath = "./resources/human/GDC/TCGA/COAD/"  
project.resultsPath = "./results2/human/GDC/TCGA/"
#*******************************************************************************************************
project.settings.DEA.H0.Log2FC = 1
project.settings.DEA.padj = 0.05
project.settings.DEA.cutoff.log2FC = 1
#*******************************************************************************************************
project.rnaseq.runHTSEQ = FALSE
prject.rnaseq.runKallisto = FALSE
project.rnaseq.runDEXSeq = FALSE
project.rnaseq.run.da = FALSE
project.rnaseq.run.go = FALSE
project.rnaseq.run.expressionSurvival = FALSE
project.rnaseq.run.immune = TRUE
project.rnaseq.run.immune.survival = FALSE
#*******************************************************************************************************
prject.rnaseq.ensemblSpecie = "hsapiens_gene_ensembl" 
project.formula = ~GROUP
project.groupsToCompare = list(c("GROUP","Primary Tumor","Solid Tissue Normal"))  #Treatment vs control

project.groupsToCompare.inmune = list(c("Primary Tumor","Solid Tissue Normal"))  #Treatment vs control
project.groupsToCompare.inmune.orderFactor = c("Solid Tissue Normal","Primary Tumor")
project.groupsToCompare.inmune.size.cibersort = c(width=15,height = 20)
project.groupsToCompare.inmune.colCount.cibersort = 4

project.groupsToCompare.inmune.size.xcell = c(width=10,height = 65)
project.groupsToCompare.inmune.colCount.xcell = 4


project.groupsToCompare.inmune.size.quantiseq = c(width=10,height = 15)
project.groupsToCompare.inmune.colCount.quantiseq = 4

project.groupsToCompare.inmune.size.epic = c(width=10,height = 12)
project.groupsToCompare.inmune.colCount.epic = 4
#*******************************************************************************************************
geneListFilterPath = "./resources/GlycoGeneList/Glyco_HS_GeneList.csv"
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList_filtered.csv"
samplesListFPKMPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSeq-FPKM_samples.csv"
#*******************************************************************************************************
samplesList = HELPER_LoadCSV(samplesListPath)
samplesList = samplesList[!(is.na(samplesList$GROUP) | samplesList$GROUP == "not reported"), ]
levels(samplesList$GROUP)
samplesList$GROUP =  droplevels(samplesList$GROUP)
levels(samplesList$GROUP)
project.sample.name = samplesList$Sample_ID
project.sample.group = as.factor(samplesList$GROUP)
table(samplesList$GROUP)


source("./src/shared/settings.r")

