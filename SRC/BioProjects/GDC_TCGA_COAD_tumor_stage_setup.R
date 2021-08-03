source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "TCGA"
project.specie ="HUMAN"
project.name = "COAD_TUMOR_STAGE"   
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
project.formula = ~ tumor_stage
project.groupsToCompare = list(c("tumor_stage","stage iv","stage iii"),
                               c("tumor_stage","stage iv","stage ii"),
                               c("tumor_stage","stage iv","stage i"),
                               c("tumor_stage","stage iii","stage ii"),
                               c("tumor_stage","stage iii","stage i"),
                               c("tumor_stage","stage ii","stage i"),
                               c("tumor_stage","stage iv","Solid Tissue Normal"),
                               c("tumor_stage","stage iii","Solid Tissue Normal"),
                               c("tumor_stage","stage ii","Solid Tissue Normal"),
                               c("tumor_stage","stage i","Solid Tissue Normal"))  #Treatment vs control

project.groupsToCompare.inmune = list(c("stage iv","stage iii"),
                                      c("stage iv","stage ii"),
                                      c("stage iv","stage i"),
                                      c("stage iii","stage ii"),
                                      c("stage iii","stage i"),
                                      c("stage ii","stage i"),
                                      c("stage iv","Solid Tissue Normal"),
                                      c("stage iii","Solid Tissue Normal"),
                                      c("stage ii","Solid Tissue Normal"),
                                      c("stage i","Solid Tissue Normal"))  #Treatment vs control
project.groupsToCompare.inmune.orderFactor = c("Solid Tissue Normal","stage i","stage ii","stage iii","stage iv")

project.groupsToCompare.inmune.size.cibersort = c(width=15,height = 20)
project.groupsToCompare.inmune.colCount.cibersort = 4

project.groupsToCompare.inmune.size.xcell = c(width=15,height = 90)
project.groupsToCompare.inmune.colCount.xcell = 4

project.groupsToCompare.inmune.size.quantisec = c(width=15,height = 30)
project.groupsToCompare.inmune.colCount.quantiseq = 4

project.groupsToCompare.inmune.size.epic = c(width=15,height = 25)
project.groupsToCompare.inmune.colCount.epic = 4

#*******************************************************************************************************
geneListFilterPath = "./resources/GlycoGeneList/Glyco_HS_GeneList.csv"
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList_filtered.csv"
samplesListFPKMPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSeq-FPKM_samples.csv"
#*******************************************************************************************************
samplesList = HELPER_LoadCSV(samplesListPath)
samplesList$tumor_stage = as.character(samplesList$tumor_stage)

samplesList$tumor_stage = as.factor(samplesList$tumor_stage)
levels(samplesList$tumor_stage)
samplesList$tumor_stage =  droplevels(samplesList$tumor_stage)
levels(samplesList$tumor_stage)
project.sample.name = samplesList$Sample_ID
project.sample.group = as.factor(samplesList$tumor_stage)

table(project.sample.group)

#*******************************************************************************************************
source("./src/shared/settings.r")





