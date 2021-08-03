#https://wiki.nci.nih.gov/display/TCGA/Clinical+Data+Overview

source("./src/shared/packages.r")
#************************************************************************************************************************************************************
project.type = "Clinical"
project.resultsPath = "./results/human/GDC/TCGA/Clinical/"
clinicalPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_ClinicalMerge.csv"
samplesListPath = "./resources/human/GDC/TCGA/COAD/TCGA_COAD_HTSEQ_Clinical_SampleList_filtered.csv"

outputDir = project.resultsPath
dir.create(file.path(outputDir), recursive = TRUE, showWarnings = FALSE)
source("./src/shared/settings.r")
#********************************************#****************************************************************************************************************
# Variable descriptive analysis
#****************************************************

clinicalData = HELPER_LoadCSV(clinicalPath)
samplesList = HELPER_LoadCSV(samplesListPath)
#filter by selected samples ID
clinicalData  = clinicalData[clinicalData$submitter_id %in% samplesList$Case_ID,]

HELPER_SAVE_DATA_FRAME(clinicalData, paste0(outputDir,"Clinical.csv")) 
#TODO esto ver del mail que hay una forma más limpia
clinicalData[clinicalData=="--" | clinicalData == "not reported"] = NA
clinicalData = droplevels(clinicalData)
clinicalData$age_at_diagnosis = as.numeric(clinicalData$age_at_diagnosis)/365
clinicalData$tumor_stage = gsub("ia","i",clinicalData$tumor_stage)
clinicalData$tumor_stage = gsub("ib","i",clinicalData$tumor_stage)
clinicalData$tumor_stage = gsub("ic","i",clinicalData$tumor_stage)
clinicalData$tumor_stage = gsub("va","v",clinicalData$tumor_stage)
clinicalData$tumor_stage = gsub("vb","v",clinicalData$tumor_stage)

clinicalData$days_to_birth = abs(clinicalData$days_to_birth/365)
colnames(clinicalData)[21] = "Age"
clinicalData = clinicalData[!is.na(clinicalData$vital_status),]
#*************
#missing values
#*************
a = aggr(clinicalData, plot = FALSE)
b = as.data.frame(a$missings)

pdf(paste0(outputDir,"missingVars.pdf"))
barplot(b$Count, names.arg = b$Variable, las=2, col = "red")
dev.off()
#*************
#plots
#*************

tumorStage= data.frame(clinicalData$gender,clinicalData$tumor_stage)
barplotGraphic(tumorStage, beside = TRUE,title ="Tumor Stage", yTitle = "Frecuency", path = paste0(outputDir,"tumorStageByGender.pdf"))


diagnosisAgeByGender = data.frame(GENDER = clinicalData$gender, DIAGNOSIS_AGE = clinicalData$age_at_diagnosis, stringsAsFactors = FALSE)
barplotGraphic(diagnosisAgeByGender,breaks= seq(30,100,10),  beside = TRUE,title ="Diagnosis Age", yTitle = "Frecuency",xTitle = "Years",path = paste0(outputDir,"diagnosisAgeByGender.pdf"))

bmiByGender = data.frame(GENDER = clinicalData$gender, BMI = clinicalData$bmi, stringsAsFactors = FALSE)
barplotGraphic(bmiByGender,labels= c("Underweight","Normal","Overweight","Obesity"),breaks= c(0,18.5,24.9,29.9,40),  beside = TRUE,title ="BMI", yTitle = "Frecuency",path = paste0(outputDir,"BMIByGender.pdf"))


weightByGender = data.frame(GENDER = clinicalData$gender, WEIGHT = clinicalData$weight, stringsAsFactors = FALSE)
barplotGraphic(weightByGender,breaks =  c(30,40,50,60,70,80,90,100,110,200), beside = TRUE,title ="Weigth", yTitle = "Frecuency",path = paste0(outputDir,"weigthByGender.pdf"))

pieGraphic(clinicalData$gender,title = "Gender",path = paste0(outputDir,"gender.pdf"))
pieGraphic(clinicalData$race,title = "Race",path = paste0(outputDir,"race.pdf"))
pieGraphic(clinicalData$morphology,title = "Morphology",path = paste0(outputDir,"morphology.pdf"))
pieGraphic(clinicalData$tissue_or_organ_of_origin,title = "Tissue or organ of origin",path = paste0(outputDir,"tissue.pdf"))
pieGraphic(clinicalData$tumor_stage[clinicalData$vital_status == "dead"],title = "Dead - Tumor stage",path = paste0(outputDir,"Dead-Stage.pdf"))
pieGraphic(clinicalData$tumor_stage[clinicalData$vital_status == "alive"],title = "Alive - Tumor stage",path = paste0(outputDir,"Alive-Stage.pdf"))

#************************************************************************************************************************************************************************
#SURVIVAL
#************************************************************************************************************************************************************************
survivalData = GET_SURVIVAL_TIME(timeToLastFollow = clinicalData$days_to_last_follow_up,
                                 vitalStatus = clinicalData$vital_status, 
                                 daysToDeath = clinicalData$days_to_death,
                                 otherVars = clinicalData)
survivalData = survivalData[,-c(3,4,5,16,18,19,20,21)]

GET_OVERALL_SURVIVAL(survivalData = survivalData,
                     outputDir = paste0(outputDir,"survival/"))


GET_KAPLAN_MEIER_PLOT(survivalData = survivalData,
              outputDir = paste0(outputDir, "survival/Kaplan-Meier/"))

survivalData = survivalData[,c(1,2,3,4,11,12,15,30,40,46,47,53,54)]
survivalData = na.omit(survivalData)
GET_COXPH_PLOT(survivalData = survivalData,
               outputDir = paste0(outputDir, "survival/COXPH/"))



