source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "TCGA"
project.specie ="HUMAN"
project.name = "COAD_LOST_REPAIR_PROTEINS_EXPRESSION"   
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
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList.csv"
samplesListFPKMPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSeq-FPKM_samples.csv"
#*******************************************************************************************************
samplesList = HELPER_LoadCSV(samplesListPath)
samplesList = samplesList[!(is.na(samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc_result) |
                is.na(samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc_result.1) | 
                is.na(samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc_result.2) |
                is.na(samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc_result.3)),]

LostRepairProtein = paste(MSH2 = samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc_result, PMSH2 = samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc_result.1, MSH6 =  samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc_result.2, MLH1 = samplesList$loss_expression_of_mismatch_repair_proteins_by_ihc_result.3, setp = "_")
LostRepairProtein = gsub(" ", "_", LostRepairProtein)
LostRepairProtein = gsub("__", "", LostRepairProtein)
LostRepairProtein = gsub("-", "_", LostRepairProtein)
LostRepairProtein[samplesList$GROUP == "Solid Tissue Normal"] = "Control"

samplesList$LostRepairProteins =  as.factor(LostRepairProtein)
levels(samplesList$LostRepairProteins)

project.formula = ~ LostRepairProteins
project.groupsToCompare = list(c("LostRepairProteins","msh2_expressed_pms2_not_expressed_msh6_expressed_mlh1_expressed","Control"),
                               c("LostRepairProteins","msh2_expressed_pms2_not_expressed_msh6_not_expressed_mlh1_expressed","Control"),
                               c("LostRepairProteins","msh2_not_expressed_pms2_not_expressed_msh6_not_expressed_mlh1_not_expressed","Control"),
                               c("LostRepairProteins","msh2_expressed_pms2_expressed_msh6_expressed_mlh1_expressed","Control"),
                              c("LostRepairProteins","msh2_expressed_pms2_not_expressed_msh6_expressed_mlh1_not_expressed","Control"),
                              c("LostRepairProteins","msh2_not_expressed_pms2_expressed_msh6_not_expressed_mlh1_expressed","Control"))




project.sample.name = samplesList$Sample_ID
project.sample.group = as.factor(samplesList$LostRepairProteins)

table(project.sample.group)

#*******************************************************************************************************
source("./src/shared/settings.r")



a = 
dim(a)
View(a)

