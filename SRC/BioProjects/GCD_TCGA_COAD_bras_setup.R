source("./src/shared/packages.r")
#*******************************************************************************************************
#Project Settings  - Setup this FIRST!
#*******************************************************************************************************
project.type = "TCGA"
project.specie ="HUMAN"
project.name = "COAD_BRAF"   
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
project.formula = ~ braf_gene_analysis_result
project.groupsToCompare = list(c("braf_gene_analysis_result","abnormal","normal"),
                               c("braf_gene_analysis_result","abnormal","Solid Tissue Normal"),
                               c("braf_gene_analysis_result","normal","Solid Tissue Normal"))
#*******************************************************************************************************
geneListFilterPath = "./resources/GlycoGeneList/Glyco_HS_GeneList.csv"
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList.csv"
samplesListFPKMPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSeq-FPKM_samples.csv"
#*******************************************************************************************************
samplesList = HELPER_LoadCSV(samplesListPath)
samplesList = samplesList[!(is.na(samplesList$braf_gene_analysis_result) | samplesList$braf_gene_analysis_result == "indeterminate"), ]
levels(samplesList$braf_gene_analysis_result)
samplesList$msi_status =  droplevels(samplesList$braf_gene_analysis_result)
levels(samplesList$braf_gene_analysis_result)
project.sample.name = samplesList$Sample_ID
project.sample.group = as.factor(samplesList$braf_gene_analysis_result)

table(project.sample.group)

#*******************************************************************************************************
source("./src/shared/settings.r")





