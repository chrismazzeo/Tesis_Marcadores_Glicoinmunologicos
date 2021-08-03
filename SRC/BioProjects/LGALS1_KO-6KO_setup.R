source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "RNASEQ"
project.specie ="MOUSE"
project.name = "LGALS1_KO-6KO_v2"   
project.resourcePath = "./resources/mouse/LGALS1_KO-6KO/"  
project.resultsPath = "./results/mouse/"
#*******************************************************************************************************
project.settings.DEA.H0.Log2FC = c(0,1)
project.settings.DEA.padj = 0.05
project.settings.DEA.cutoff.log2FC = 1
#*******************************************************************************************************
project.rnaseq.runHTSEQ = TRUE
prject.rnaseq.runKallisto = FALSE
project.rnaseq.runDEXSeq = FALSE
project.rnaseq.run.da = FALSE
project.rnaseq.run.go = FALSE
project.rnaseq.run.expressionSurvival = FALSE
project.rnaseq.run.immune = FALSE
project.rnaseq.run.immune.survival = FALSE
#*******************************************************************************************************
prject.rnaseq.ensemblSpecie = "mmusculus_gene_ensembl"
#*******************************************************************************************************
project.immucc.svrPath =  "./resources/mouse/LGALS1_KO-6KO/immucc/LGALS1_KO-6KO.SVR.csv"   
project.immucc.llsrPath = "./resources/mouse/LGALS1_KO-6KO/immucc/LGALS1_KO-6KO.LLSR.csv"   
#*******************************************************************************************************
geneListFilterPath = "./resources/GlycoGeneList/Glyco_MM_GeneList.csv"
samplesListPath = "./resources/mouse/LGALS1_KO-6KO/LGALS1_KO-6KO_sampleList.csv"
#*******************************************************************************************************
samplesList = HELPER_LoadCSV(samplesListPath)
project.formula = ~ GROUP
project.sample.name = samplesList$SAMPLE
project.sample.group = as.factor(samplesList$GROUP)
project.groupsToCompare = list(c("GROUP","LGALS1_KO","AOM_DSS"))  #Treatment vs control

# samplesList = samplesList[samplesList$GROUP == "LGALS1_KO",]
# project.formula = ~ TUMOR_CLASS
# project.sample.name = samplesList$SAMPLE
# project.sample.group = as.factor(samplesList$TUMOR_CLASS)
# project.groupsToCompare = list(c("TUMOR_CLASS","MED","LOW"))  #Treatment vs control
#*******************************************************************************************************
source("./src/shared/settings.r")



