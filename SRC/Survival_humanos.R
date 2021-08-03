#************************************************************************************************************************************************************************
#SURVIVAL
#************************************************************************************************************************************************************************
source("./src/shared/packages.r")
library(readr)
library(survival)
library(survminer)

survivalData <- read_delim("/Volumes/Externo/Google Drive/Bioinformática/PFI/Tesina/Tesis Final/Resultados Finales/Anexos/Anexo Survival/Anexos Tesis - Survival.csv", 
                          ";", escape_double = FALSE, trim_ws = TRUE)


outputDir = "./../"

#nos quedamos solo con los tumores
survivalData = survivalData[survivalData$Grupo=="Tumor",]

#nota el vital status tiene que ser alive o dead
survivalData = GET_SURVIVAL_TIME(timeToLastFollow = survivalData$`Días hasta el último seguimiento`,
                                 vitalStatus = survivalData$Estado, 
                                 daysToDeath = survivalData$`Días hasta la muerte`,
                                 otherVars = survivalData)

#limpiamos las columnas que no necesitamos mas
survivalData = survivalData[,-c(3,4,5,9,10,11)]

#sacamos white spaces y otros caracteres que rompen
colnames(survivalData) = gsub(" ", "_", colnames(survivalData))
colnames(survivalData) = gsub("\\+", "", colnames(survivalData))
colnames(survivalData) = gsub("-", "_", colnames(survivalData))

survivalData$Time = as.numeric(survivalData$Time)
View(survivalData)


#overall survival
GET_OVERALL_SURVIVAL(survivalData = survivalData,
                     outputDir = paste0(outputDir,"survival/"))

#kaplan meier de cada variable
GET_KAPLAN_MEIER_PLOT(survivalData = survivalData,
                      outputDir = paste0(outputDir, "survival/Kaplan-Meier/"))
#************************************************************
#hazard ratio
#************************************************************
# Si una categoria no tiene eventos va a romper, por eso tenemos que ir viendo que variables son las que tienen alguna categoria sin eventos
# para continuas no se porque pero a veces tambien rompe
#columnas que rompen siempre 5, 22, 29, 30,32,38,41,
survivalData2 = survivalData[,c(1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,31,33,34,35,36,37,39,40)]
#survivalData2 = survivalData[,c(1,2,6)]

survivalData2 = na.omit(survivalData2)

GET_COXPH_PLOT(survivalData = survivalData2,
               outputDir = paste0(outputDir, "survival/COXPH/"))

#************************************************************
#juntamos los valores de KLB, Notum y PTN
#************************************************************
res.cut <- surv_cutpoint(survivalData, time = "Time", event = "VitalStatus",
                         variables = c(colnames(survivalData)[6:8]))

categories = surv_categorize(res.cut)[3:5]
categories
survivalData2 = data.frame(value = survivalData[,1:2], categories)
colnames(survivalData2) = c(colnames(survivalData)[1:2],colnames(categories))

all =  as.factor(paste0("Ptn: ",as.character(survivalData2$Ptn), "  Klb: ", as.character(survivalData2$Klb),"  Notum: ", as.character(survivalData2$Notum)))
survivalData2$all = all
survivalData2 = survivalData2[,c(1,2,6)]

GET_KAPLAN_MEIER_PLOT(survivalData = survivalData2,
                      outputDir = paste0(outputDir, "survival/Kaplan-Meier/"))

#solo la que tiene mas y menos supervivencia
caca = survivalData2[(survivalData2$all == "Ptn: low  Klb: low  Notum: low" | survivalData2$all == "Ptn: low  Klb: low  Notum: high"),]
caca = droplevels(caca)
GET_KAPLAN_MEIER_PLOT(survivalData = caca,
                      outputDir = paste0(outputDir, "survival/depurado/"))

#************************************************************

#Agrupado por combinaciones

 surObj = Surv(survivalData$Time, survivalData$VitalStatus)

 fit.surv <- survfit(surObj ~  Estadío  + MSI + dMMR , data = survivalData )
# 
# 
# 
# 
 ggsurv = ggsurvplot(fit.surv,   pval = TRUE, conf.in = FALSE,
          
                     surv.plot.height = 0.5, 
                     risk.table = "percentage", risk.table.y.text = FALSE, risk.table.col="strata", 
                     legend = "right",
                     surv.median.line = "hv",
                     fun =  "pct",
                     title = "",
                     xlab = "Tiempo(días)",
                     ylab = "P. Supervivencia",
                     font.legend = c(12,"plain", "black"),
       
                     tables.height = 0.25,
                     tables.theme =settings.graphics.theme,
                     cumcensor = TRUE,
                     ggtheme = theme_bw())
 ggsurv$plot
 curv_facet <- ggsurv$plot + facet_grid(Estadío ~ MSI )
 curv_facet

 



                                                                                        