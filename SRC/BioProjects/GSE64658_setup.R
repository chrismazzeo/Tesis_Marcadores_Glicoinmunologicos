source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "MICROARRAY"
project.specie ="MOUSE"
project.name = "GSE64658"   
project.resourcePath = "./resources/mouse/"  
project.resultsPath = "./results/mouse/"
#*******************************************************************************************************
project.microarray.cel.source = "GEO"  #ExpressArray  #Local
project.microarray.downloadAccessionNumber = "GSE64658"
project.microarray.backgroundCorrectionMode = "GCRMA"  #"RMA"
project.microarray.annotationDB = mouse4302.db  #annotation(affyData)  #esto tengo que ver como como lo asigno automaticamente a la variable
project.microarray.settings.lowIntensityFilter = 2.5
project.settings.DEA.H0.Log2FC = 0
project.settings.DEA.padj = 0.05
project.settings.DEA.cutoff.log2FC = 1
#*******************************************************************************************************
project.immucc.svrPath =  "./resources/mouse/GSE64658/immucc/GSE64658.SVR.csv"   #Affy 4302 
project.immucc.llsrPath = "./resources/mouse/GSE64658/immucc/GSE64658.LLSR.csv"   #Affy 4302 
#*******************************************************************************************************
project.groupsToCompare =  c("Distal_Colon_AOM_DSS_Early-Distal_Colon_Control",
                             "Distal_Colon_AOM_DSS_Late-Distal_Colon_Control",
                             "Proximal_Colon_AOM_DSS_Early-Proximal_Colon_Control",
                             "Proximal_Colon_AOM_DSS_Late-Proximal_Colon_Control",
                             "Distal_Colon_AOM_DSS_Late-Distal_Colon_AOM_DSS_Early",
                             "Proximal_Colon_AOM_DSS_Late-Proximal_Colon_AOM_DSS_Early")

geneListFilterPath = "./resources/GlycoGeneList/Glyco_MM_GeneList.csv"
samplesListPath = "./resources/mouse/GSE64658/GSE64658_SampleList.csv"
#*******************************************************************************************************
source("./src/shared/settings.r")

