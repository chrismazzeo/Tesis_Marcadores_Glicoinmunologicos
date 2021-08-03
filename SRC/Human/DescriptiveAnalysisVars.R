source("./src/functions/unzipHTCounts.r")
source("./src/functions/graphics.r")


list.of.packages <- c("psych","pastecs","magrittr","kableExtra","survival","plotly","readxl","XML","dplyr","stringr","survminer","VIM")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(psych)
library(pastecs)
library(knitr)
#install.packages("magrittr")
library(magrittr)
#install.packages("kableExtra")
library(kableExtra)
library(survival)
library(readr)
library(plotly)
library(readxl)
library(XML)
library(dplyr)
library(stringr)
library(survminer)
library(VIM)
#*******


#*************
#Load clinical
#*************
clinical = loadTSV(output.clinicalSamples)
clinical[clinical=="--"] = NA
clinical$age_at_diagnosis = as.numeric(clinical$age_at_diagnosis)/365
clinical$tumor_stage = gsub("ia","i",clinical$tumor_stage)
clinical$tumor_stage = gsub("ib","i",clinical$tumor_stage)
clinical$tumor_stage = gsub("ic","i",clinical$tumor_stage)
clinical$tumor_stage = gsub("va","v",clinical$tumor_stage)
clinical$tumor_stage = gsub("vb","v",clinical$tumor_stage)
#*************
#missing values
#*************
a =aggr(clinical, plot = FALSE)
b = as.data.frame(a$missings)

pdf(paste0(folder.output.da,"missingVars.pdf"))
barplot(b$Count, names.arg = b$Variable, las=2, col = "red")
dev.off()
#*************
#plots
#*************

tumorStage= data.frame(clinical$gender,clinical$tumor_stage)
barplotGraphic(tumorStage, beside = TRUE,title ="Tumor Stage", yTitle = "Frecuency", path = paste0(folder.output.da,"tumorStageByGender.pdf"))


diagnosisAgeByGender = data.frame(GENDER = clinical$gender, DIAGNOSIS_AGE = clinical$age_at_diagnosis, stringsAsFactors = FALSE)
barplotGraphic(diagnosisAgeByGender,breaks= seq(30,100,10),  beside = TRUE,title ="Diagnosis Age", yTitle = "Frecuency",xTitle = "Years",path = paste0(folder.output.da,"diagnosisAgeByGender.pdf"))

bmiByGender = data.frame(GENDER = clinical$gender, BMI = clinical$bmi, stringsAsFactors = FALSE)
barplotGraphic(bmiByGender,labels= c("Underweight","Normal","Overweight","Obesity"),breaks= c(0,18.5,24.9,29.9,40),  beside = TRUE,title ="BMI", yTitle = "Frecuency",path = paste0(folder.output.da,"BMIByGender.pdf"))


weightByGender = data.frame(GENDER = clinical$gender, WEIGHT = clinical$weight, stringsAsFactors = FALSE)
barplotGraphic(weightByGender,breaks =  c(30,40,50,60,70,80,90,100,110,200), beside = TRUE,title ="Weigth", yTitle = "Frecuency",path = paste0(folder.output.da,"weigthByGender.pdf"))

pieGraphic(clinical$gender,title = "Gender",path = paste0(folder.output.da,"gender.pdf"))
pieGraphic(clinical$race,title = "Race",path = paste0(folder.output.da,"race.pdf"))
pieGraphic(clinical$morphology,title = "Morphology",path = paste0(folder.output.da,"morphology.pdf"))
pieGraphic(clinical$tissue_or_organ_of_origin,title = "Tissue or organ of origin",path = paste0(folder.output.da,"tissue.pdf"))
pieGraphic(clinical$tumor_stage[clinical$vital_status == "dead"],title = "Dead - Tumor stage",path = paste0(folder.output.da,"Dead-Stage.pdf"))
pieGraphic(clinical$tumor_stage[clinical$vital_status == "alive"],title = "Alive - Tumor stage",path = paste0(folder.output.da,"Alive-Stage.pdf"))
pieGraphic(clinical$microsatellite_instability,title = "Microsatellite inestability",path = paste0(folder.output.da,"microsatelliteInestability.pdf"))
pieGraphic(clinical$mononucleotide_and_dinucleotide_marker_panel_analysis_status,title = "MSI status",path = paste0(folder.output.da,"MSIStatus.pdf"))
pieGraphic(clinical$braf_gene_analysis_result, title = "Braf Analysis",path = paste0(folder.output.da,"BrafAnalysis.pdf")) 
pieGraphic(clinical$kras_mutation_found, title = "Kras Analysis",path = paste0(folder.output.da,"KrasAnalysis.pdf")) 
pieGraphic(clinical$loss_expression_of_mismatch_repair_proteins_by_ihc, title = "Loss expression of mismatch repair protein",path = paste0(folder.output.da,"IHCLoss.pdf")) 
pieGraphic(clinical$loss_expression_of_mismatch_repair_proteins_by_ihc_result, title = "Loss expression of mismatch repair protein by Ihc",path = paste0(folder.output.da,"IHC1.pdf")) 
pieGraphic(clinical$.loss_expression_of_mismatch_repair_proteins_by_ihc_result.2., title = "Loss expression of mismatch repair protein by Ihc",path = paste0(folder.output.da,"IHC2.pdf")) 
pieGraphic(clinical$.loss_expression_of_mismatch_repair_proteins_by_ihc_result.3., title = "Loss expression of mismatch repair protein by Ihc",path = paste0(folder.output.da,"IHC3.pdf")) 
pieGraphic(clinical$.loss_expression_of_mismatch_repair_proteins_by_ihc_result.4., title = "Loss expression of mismatch repair protein by Ihc",path = paste0(folder.output.da,"IHC4.pdf")) 

#******************************************************************************************************************************
#Survival   kaplan-meier
#******************************************************************************************************************************
t = data.frame(Time = as.numeric(clinical$days_to_last_follow_up),Status = clinical$vital_status)
t$Time[!is.na(clinical$days_to_death)] = as.numeric(clinical$days_to_death[!is.na(clinical$days_to_death)])
index = t$Status == "alive"
t$Status = 1
t$Status[index] = 0
clinical = data.frame(clinical, Time = t$Time, Status = t$Status)

#removemos los que tengan tiempo 0
clinical = clinical[clinical$Time>0,]
dim(clinical)



pdf(paste0(folder.output.da,"survival.pdf"))
t = data.frame (Time = clinical$Time, Status = clinical$Status, Stage = clinical$tumor_stage, MSI = clinical$mononucleotide_and_dinucleotide_marker_panel_analysis_status, stringsAsFactors = TRUE)
t = na.omit(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ 1, data = t )
ggsurvplot(fit1, data = t,  conf.in =TRUE, pval = TRUE, risk.table = "abs_pct" )



fit1 <- survfit(surObj ~Stage , data = t )
ggsurvplot(fit1, data = t,  pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  Stage , data = t)
ggforest(fit.coxph, data = t)

fit1 <- survfit(surObj ~ MSI, data = t )
ggsurvplot(fit1, data = t,   pval = TRUE )
fit.coxph <- coxph(surObj  ~  MSI , data = t)
ggforest(fit.coxph, data = t)

t = data.frame (Time = clinical$Time, Status = clinical$Status,KRAS_Mutation = clinical$kras_mutation_found, stringsAsFactors = TRUE)
t = na.omit(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ KRAS_Mutation, data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  KRAS_Mutation , data = t)
ggforest(fit.coxph, data = t)


t = data.frame (Time = clinical$Time, Status = clinical$Status,BRAF_Mutation = clinical$braf_gene_analysis_result, stringsAsFactors = TRUE)
t = na.omit(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ BRAF_Mutation, data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  BRAF_Mutation , data = t)
ggforest(fit.coxph, data = t)



t = data.frame (Time = clinical$Time, Status = clinical$Status,MLH1 = clinical$loss_expression_of_mismatch_repair_proteins_by_ihc_result, stringsAsFactors = TRUE)
t = na.omit(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ MLH1, data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  MLH1 , data = t)
ggforest(fit.coxph, data = t)


t = data.frame (Time = clinical$Time, Status = clinical$Status,MLH2 = clinical$.loss_expression_of_mismatch_repair_proteins_by_ihc_result.2., stringsAsFactors = TRUE)
t = na.omit(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ MLH2, data = t )
ggsurvplot(fit1, data = t,  pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  MLH2 , data = t)
ggforest(fit.coxph, data = t)


t = data.frame (Time = clinical$Time, Status = clinical$Status,PMS2 = clinical$.loss_expression_of_mismatch_repair_proteins_by_ihc_result.3., stringsAsFactors = TRUE)
t = na.omit(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ PMS2, data = t )
ggsurvplot(fit1, data = t,  pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  PMS2 , data = t)
ggforest(fit.coxph, data = t)


t = data.frame (Time = clinical$Time, Status = clinical$Status,MSH6 = clinical$.loss_expression_of_mismatch_repair_proteins_by_ihc_result.4., stringsAsFactors = TRUE)
t = na.omit(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ MSH6, data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  MSH6 , data = t)
ggforest(fit.coxph, data = t)


t = data.frame (Time = clinical$Time, Status = clinical$Status,Lyn = clinical$lymphatic_invasion, stringsAsFactors = TRUE)
t = na.omit(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ Lyn, data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  Lyn , data = t)
ggforest(fit.coxph, data = t)


t = data.frame (Time = clinical$Time, Status = clinical$Status,Rx = clinical$radiation_therapy, stringsAsFactors = TRUE)
t = na.omit(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ Rx, data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  Rx , data = t)
ggforest(fit.coxph, data = t)


t = data.frame (Time = clinical$Time, Status = clinical$Status,Lyn = clinical$lymphatic_invasion,stringsAsFactors = TRUE)
t = na.omit(t)
dim(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ LynpathicInvasion, data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  Lyn , data = t)
ggforest(fit.coxph, data = t)

t = data.frame (Time = clinical$Time, Status = clinical$Status,RxTherapy = clinical$radiation_therapy,stringsAsFactors = TRUE)
t = na.omit(t)
dim(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ RxTherapy, data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )


t = data.frame (Time = clinical$Time, Status = clinical$Status,Polyps = clinical$colon_polyps_present,stringsAsFactors = TRUE)
t = na.omit(t)
dim(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ Polyps, data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~  Polyps , data = t)
ggforest(fit.coxph, data = t)


t = data.frame (Time = clinical$Time, Status = clinical$Status,Lyn = clinical$lymphatic_invasion, Rx = clinical$radiation_therapy,stringsAsFactors = TRUE)
t = na.omit(t)
dim(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ Lyn+ Rx , data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~   Lyn+ Rx , data = t)
ggforest(fit.coxph, data = t)


t = data.frame (Time = clinical$Time, Status = clinical$Status,MLH1 = clinical$loss_expression_of_mismatch_repair_proteins_by_ihc_result,MLH2 = clinical$.loss_expression_of_mismatch_repair_proteins_by_ihc_result.2.,PMS2 = clinical$.loss_expression_of_mismatch_repair_proteins_by_ihc_result.3., MSH6 = clinical$.loss_expression_of_mismatch_repair_proteins_by_ihc_result.4.,stringsAsFactors = TRUE)
t = na.omit(t)
dim(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ MLH1+ MLH2 + PMS2 + MSH6 , data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )

fit.coxph <- coxph(surObj  ~   MLH1+ MLH2 + PMS2 + MSH6 , data = t)
ggforest(fit.coxph, data = t)

t = data.frame (Time = clinical$Time, Status = clinical$Status,VenousInvasion = clinical$venous_invasion,stringsAsFactors = TRUE)
t = na.omit(t)
dim(t)
surObj = Surv(time = t$Time, event = t$Status)
fit1 <- survfit(surObj ~ VenousInvasion , data = t )
ggsurvplot(fit1, data = t, pval = TRUE, risk.table = "abs_pct" )
fit.coxph <- coxph(surObj  ~   VenousInvasion , data = t)
ggforest(fit.coxph, data = t)

#nota: si meto todas las variables juntas hay demsiados NA y no tiene potencia
dev.off()








