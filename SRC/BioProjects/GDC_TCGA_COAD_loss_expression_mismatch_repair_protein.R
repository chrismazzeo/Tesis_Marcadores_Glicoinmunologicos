source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "TCGA"
project.specie ="HUMAN"
project.name = "COAD_Loss_expression_mismatch_repair_proteins"   
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
project.rnaseq.run.immune.survival = FALSE
#*******************************************************************************************************
prject.rnaseq.ensemblSpecie = "hsapiens_gene_ensembl" 
project.formula = ~ loss_expression_of_mismatch_repair_proteins_by_ihc
project.groupsToCompare = list(c("loss_expression_of_mismatch_repair_proteins_by_ihc","yes","no"),
                               c("loss_expression_of_mismatch_repair_proteins_by_ihc","yes","Solid Tissue Normal"),
                               c("loss_expression_of_mismatch_repair_proteins_by_ihc","no","Solid Tissue Normal"))

project.groupsToCompare.inmune = list(c("yes","no"),
                                      c("yes","Solid Tissue Normal"),
                                      c("no","Solid Tissue Normal"))  #Treatment vs control
project.groupsToCompare.inmune.orderFactor = c("Solid Tissue Normal","no","yes")
project.groupsToCompare.inmune.size.cibersort = c(width=10,height = 23)
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

levels(samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc)
samplesList$msi_status =  droplevels(samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc)
levels(samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc)
project.sample.name = samplesList$Sample_ID
project.sample.group = as.factor(samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc)

table(project.sample.group)

#*******************************************************************************************************
source("./src/shared/settings.r")





