source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "RNASEQ"
project.specie ="MOUSE"
project.name = "GSE57533"   
project.resourcePath = "./resources/mouse/GSE57533/"  
project.resultsPath = "./results/mouse/"
#*******************************************************************************************************
project.immucc.svrPath = NA
project.immucc.llsrPath = NA 
project.settings.DEA.H0.Log2FC = 0
project.settings.DEA.padj = 0.05
project.settings.DEA.cutoff.log2FC = 1
#*******************************************************************************************************
project.immucc.svrPath =  "./resources/mouse/GSE57533/immucc/GSE57533.SVR.csv"   
project.immucc.llsrPath = "./resources/mouse/GSE57533/immucc/GSE57533.LLSR.csv"   
#*******************************************************************************************************
prject.rnaseq.ensemblSpecie = "mmusculus_gene_ensembl"
project.groupsToCompare = list(c("GROUP","DSS","CONTROL"),
                               c("GROUP","AOM","CONTROL"),
                               c("GROUP","AOM","DSS"))  #Treatment vs control

geneListFilterPath = "./resources/GlycoGeneList/Glyco_MM_GeneList.csv"
samplesListPath = "./resources/mouse/GSE57533/GSE57533_sampleList.csv"
#*******************************************************************************************************
project.rnaseq.runHTSEQ = TRUE
prject.rnaseq.runKallisto = TRUE
project.rnaseq.runDEXSeq = TRUE
project.rnaseq.run.immune = TRUE
#*******************************************************************************************************
source("./src/shared/settings.r")


