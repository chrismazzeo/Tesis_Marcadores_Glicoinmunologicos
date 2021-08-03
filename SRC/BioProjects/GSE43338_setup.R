source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "MICROARRAY"
project.specie ="MOUSE"
project.name = "GSE43338"   
project.resourcePath = "./resources/mouse/"  
project.resultsPath = "./results/mouse/"
#*******************************************************************************************************
project.microarray.cel.source = "GEO"  #ExpressArray  #Local
project.microarray.downloadAccessionNumber = "GSE43338"
project.microarray.backgroundCorrectionMode = "GCRMA"  #"RMA"
project.microarray.annotationDB = moe430a.db  #annotation(affyData) 
project.microarray.settings.lowIntensityFilter = 2.5
project.immucc.svrPath = NA
project.immucc.llsrPath = NA 
project.settings.DEA.H0.Log2FC = 0
project.settings.DEA.padj = 0.05
project.settings.DEA.cutoff.log2FC = 1
#*******************************************************************************************************
project.immucc.svrPath = NA
project.immucc.llsrPath = NA 
#*******************************************************************************************************
project.groupsToCompare =  c("colorectal_tumor_colitis_associated-colorectal_control_epithelium_colitis_associated",
                             "colorectal_tumor_sporadic-colorectal_control_epithelium_sporadic",
                             "colorectal_tumor_colitis_associated-colorectal_tumor_sporadic",
                             "colorectal_control_epithelium_colitis_associated-colorectal_control_epithelium_sporadic"
                             )

geneListFilterPath = "./resources/GlycoGeneList/Glyco_MM_GeneList.csv"
samplesListPath = "./resources/mouse/GSE43338/GSE43338_SampleList.csv"
#*******************************************************************************************************
source("./src/shared/settings.r")




