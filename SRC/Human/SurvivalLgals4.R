library(readr)

Samples <- read_delim("Resources/GDC/TCGA/COAD/Reads/GDC-TCGA-COAD-HTSeq-RawCounts-SAMPLES.tsv", 
                                                    "\t", escape_double = FALSE, trim_ws = TRUE)
Clinical <- read_delim("Resources/GDC/TCGA/COAD/Clinical/GDC-TCGA-COAD-CLINICAL.tsv", 
                                                            "\t", escape_double = FALSE, trim_ws = TRUE)

library(tidyverse)
library(stringr)
#seleccionamos que tipo de tejido
a = Samples [Samples$`Sample Type` == "Primary Tumor" | Samples$`Sample Type` == "Solid Tissue Normal",]
#removemos top down
a = Samples %>% filter(str_detect(Samples$`Sample ID`,"-01A") | str_detect(Samples$`Sample ID`,"-11A"))

#removemos replicas
a= a[!duplicated(a$`Sample ID`),]



#b = merge(x = Clinical, y = Samples, by = NULL)
#b = b[b$submitter_id == b$`Case ID`,]
#b = b[!duplicated(b$`Sample ID`),]
#b = b %>% filter(str_detect(b$`Sample ID`,"-01A") | str_detect(b$`Sample ID`,"-11A"))
#View(b)
#dim(b)

CountMatrix <- read_csv("/Volumes/Externo/Projects/R pipelines/PFI ROADMAP/Results/GDC_TCGA_COAD/CountMatrix/GDC_TCGA_COAD_CountMatrix.csv")


CountMatrix = as.data.frame(GDC_TCGA_COAD_CountMatrix)
LGALS4 = CountMatrix[CountMatrix$SYMBOL %in% "LGALS4",]

LGALS4 = LGALS4[1, 4:dim(LGALS4)[2]]
dim(LGALS4)
LGALS4 = t(LGALS4)
dim(LGALS4)

colnames(LGALS4) = "LGALS4 Values"
LGALS4  = as.data.frame(LGALS4)

LGALS4 = cbind(LGALS4,Type = Samples$`Sample Type`, caseID = Samples$`Case ID`)
dim(LGALS4)

#remove duplicated
LGALS4 = LGALS4[LGALS4$Type == "Primary Tumor",]
LGALS4 = LGALS4[!duplicated(LGALS4$caseID),]f
dim(LGALS4)
dim(Clinical)
summary(LGALS4)
quantile(LGALS4$`LGALS4 Values`)


aca tengo que ver bien como mergear,tengo que armarme algo, porque estoy siempre en la misma duda
Podria levantar toda una base de datos y dejarme de joder
  

merge(x = LGALS4, y = Clinical,by.x= LGALS4$caseID, by.y = Clinical$submitter_id)

LGALS4 = cbind(LGALS4,DaysToDeath = Clinical$days_to_death, Status = Clinical$vital_status)
dim(LGALS4)
LGALS4 = LGALS4[LGALS4$Status == "alive",]
dim(LGALS4)
LGALS4 = LGALS4[LGALS4$DaysToDeath != "--",]

dim(LGALS4)

#install.packages("varhandle")
library(varhandle)
LGALS4$Days = as.numeric(unfactor(LGALS4$Days))
LGALS4 = LGALS4[LGALS4$Days != 0,]
dim(LGALS4)



library(survival) 
KM0 <- survfit(Surv(Days, `LGALS4 Values`) ~ 1,  type="kaplan-meier", conf.type="log", data=LGALS4)
plot(KM0,lty = 2:3, col = c(100,50,200), main = "Survival rate", xlab = "Days")


boxplot(LGALS4$`LGALS4 Values`~LGALS4$Type )





plot(LGALS4$`LGALS4 Values`~LGALS4$Type, type="l")
plot(log(LGALS4$`LGALS4 Values`,2), type="l")
boxplot(log(LGALS4Sorted,2))
barplot(caca)



quanTIseq_cell_fractions <- read_delim("~/Desktop/quanTIseq_cell_fractions.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)
quanTIseq_cell_fractions = as.data.frame(quanTIseq_cell_fractions)

a =as.data.frame(t(quanTIseq_cell_fractions))

a = as.matrix(a)
View(a[2:dim(a)[1],])


barplot(a[2:dim(a)[1],], col =settings.color[1:dim(a)[2]])

