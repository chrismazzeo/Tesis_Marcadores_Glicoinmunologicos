source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "TCGA"
project.specie ="HUMAN"
project.name = "COAD_KRAS"   
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
project.rnaseq.run.immune = TRUE
project.rnaseq.run.immune.survival = TRUE
#*******************************************************************************************************
prject.rnaseq.ensemblSpecie = "hsapiens_gene_ensembl" 
project.formula = ~ kras_mutation_found
project.groupsToCompare = list(c("kras_mutation_found","yes","no"),
                              c("kras_mutation_found","yes","Solid Tissue Normal"),
                              c("kras_mutation_found","no","Solid Tissue Normal"))

#*******************************************************************************************************
geneListFilterPath = "./resources/GlycoGeneList/Glyco_HS_GeneList.csv"
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList.csv"
samplesListFPKMPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSeq-FPKM_samples.csv"
#*******************************************************************************************************
samplesList = HELPER_LoadCSV(samplesListPath)
samplesList = samplesList[!(is.na(samplesList$kras_mutation_found) | samplesList$kras_mutation_found == "indeterminate"), ]
levels(samplesList$kras_mutation_found)
samplesList$msi_status =  droplevels(samplesList$kras_mutation_found)
levels(samplesList$kras_mutation_found)
project.sample.name = samplesList$Sample_ID
project.sample.group = as.factor(samplesList$kras_mutation_found)

table(project.sample.group)

#*******************************************************************************************************
source("./src/shared/settings.r")





