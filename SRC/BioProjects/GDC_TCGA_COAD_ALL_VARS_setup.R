source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "TCGA"
project.specie ="HUMAN"
project.name = "COAD_ALL_VARS"   
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

#*******************************************************************************************************
geneListFilterPath = "./resources/GlycoGeneList/Glyco_HS_GeneList.csv"
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList_filtered.csv"
samplesListFPKMPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSeq-FPKM_samples.csv"
#*******************************************************************************************************
samplesList = HELPER_LoadCSV(samplesListPath)


samplesList$all= paste(samplesList$tumor_stage,samplesList$msi_status,samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc, sep = "__")
samplesList$all[samplesList$GROUP == "Solid Tissue Normal"] = "Control"
samplesList$all = gsub(" ","",samplesList$all)
project.sample.name = samplesList$Sample_ID
project.sample.group = as.factor(samplesList$all)
project.formula =  ~all
table(project.sample.group)


a = paste("all",levels(project.sample.group), "Control", sep ="Z")

project.groupsToCompare = strsplit(a,"Z")
project.groupsToCompare = project.groupsToCompare[2:24]
#*******************************************************************************************************
source("./src/shared/settings.r")





