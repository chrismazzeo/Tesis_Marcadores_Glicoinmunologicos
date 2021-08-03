source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "RNASEQ"
project.specie ="MOUSE"
project.name = "LGALS1_KO"   
project.resourcePath = "./resources/mouse/LGALS1_KO/"  
project.resultsPath = "./results/mouse/"
#*******************************************************************************************************
project.immucc.svrPath = NA
project.immucc.llsrPath = NA 
project.settings.DEA.H0.Log2FC = 0
project.settings.DEA.padj = 0.05
project.settings.DEA.cutoff.log2FC = 1
#*******************************************************************************************************
project.immucc.svrPath =  "./resources/mouse/LGALS1_KO/immucc/LGALS1_KO.SVR.csv"   
project.immucc.llsrPath = "./resources/mouse/LGALS1_KO/immucc/LGALS1_KO.LLSR.csv"   
#*******************************************************************************************************
prject.rnaseq.ensemblSpecie = "mmusculus_gene_ensembl"
project.groupsToCompare = list(c("GROUP","LGALS1_KO","AOM_DSS"))  #Treatment vs control

geneListFilterPath = "./resources/GlycoGeneList/Glyco_MM_GeneList.csv"
samplesListPath = "./resources/mouse/LGALS1_KO/LGALS1_KO_sampleList.csv"
#*******************************************************************************************************
project.rnaseq.runHTSEQ = FALSE
prject.rnaseq.runKallisto = FALSE
project.rnaseq.runDEXSeq = TRUE
project.rnaseq.run.immune = FALSE
#*******************************************************************************************************
source("./src/shared/settings.r")



