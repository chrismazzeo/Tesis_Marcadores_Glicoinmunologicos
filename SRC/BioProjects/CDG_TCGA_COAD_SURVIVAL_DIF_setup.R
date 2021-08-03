source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "TCGA"
project.specie ="HUMAN"
project.name = "COAD_SURVIVAL_DIF"   
project.resourcePath = "./resources/human/GDC/TCGA/COAD/"  
project.resultsPath = "./results/human/GDC/TCGA/"
#*******************************************************************************************************
project.settings.DEA.H0.Log2FC = 1
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
prject.rnaseq.ensemblSpecie = "hsapiens_gene_ensembl" 
#*******************************************************************************************************
geneListFilterPath = "./resources/GlycoGeneList/Glyco_HS_GeneList.csv"
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList_filtered.csv"
samplesListFPKMPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSeq-FPKM_samples.csv"
#*******************************************************************************************************
samplesList = HELPER_LoadCSV(samplesListPath)
samplesList = samplesList[samplesList$GROUP == "Primary Tumor",]
samplesList$SurvivalDif = ""
samplesList$SurvivalDif[samplesList$vital_status == "dead" & (samplesList$days_to_death > 100 & samplesList$days_to_death <1000)] = "Dead"
samplesList$SurvivalDif[samplesList$vital_status == "alive" & as.character(samplesList$days_to_last_follow_up) > 2000] = "Alive"
samplesList = samplesList[samplesList$SurvivalDif != "",]
table(samplesList$SurvivalDif)

project.formula = ~ SurvivalDif
project.groupsToCompare = list(c("SurvivalDif","Alive","Dead"))

project.sample.name = samplesList$Sample_ID
project.sample.group = as.factor(samplesList$SurvivalDif)

table(project.sample.group)

#*******************************************************************************************************
source("./src/shared/settings.r")

