source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "TCGA"
project.specie ="HUMAN"
project.name = "COAD_TUMOR_MSI"   
project.resourcePath = "./resources/human/GDC/TCGA/COAD/"  
project.resultsPath = "./results/human/GDC/TCGA/"
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
project.rnaseq.run.immune = FALSE
project.rnaseq.run.immune.survival = FALSE
#*******************************************************************************************************
prject.rnaseq.ensemblSpecie = "hsapiens_gene_ensembl" 
project.formula = ~ msi_status
project.groupsToCompare = list(c("msi_status","msi-h","msi-l"),
                               c("msi_status","msi-h","mss"),
                               c("msi_status","msi-h","Solid Tissue Normal"),
                               c("msi_status","msi-l","mss"),
                               c("msi_status","msi-l","Solid Tissue Normal"),
                               c("msi_status","mss","Solid Tissue Normal"))  #Treatment vs control

project.groupsToCompare.inmune = list(c("msi-h","msi-l"),
                                      c("msi-h","mss"),
                                      c("msi-h","Solid Tissue Normal"),
                                      c("msi-l","mss"),
                                      c("msi-l","Solid Tissue Normal"),
                                      c("mss","Solid Tissue Normal")) 
project.groupsToCompare.inmune.orderFactor = c("Solid Tissue Normal","mss","msi-l","msi-h")

project.groupsToCompare.inmune.size.cibersort = c(width=15,height = 20)
project.groupsToCompare.inmune.colCount.cibersort = 4

project.groupsToCompare.inmune.size.xcell = c(width=10,height = 90)
project.groupsToCompare.inmune.colCount.xcell = 4

project.groupsToCompare.inmune.size.quantisec = c(width=15,height = 30)
project.groupsToCompare.inmune.colCount.quantiseq = 4

project.groupsToCompare.inmune.size.epic = c(width=15,height = 20)
project.groupsToCompare.inmune.colCount.epic = 4

#*******************************************************************************************************
geneListFilterPath = "./resources/GlycoGeneList/Glyco_HS_GeneList.csv"
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList_filtered.csv"
samplesListFPKMPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSeq-FPKM_samples.csv"
#*******************************************************************************************************
samplesList = HELPER_LoadCSV(samplesListPath)
levels(samplesList$msi_status)
samplesList$msi_status =  droplevels(samplesList$msi_status)
levels(samplesList$msi_status)
project.sample.name = samplesList$Sample_ID
project.sample.group = as.factor(samplesList$msi_status)

table(project.sample.group)

#*******************************************************************************************************
source("./src/shared/settings.r")





