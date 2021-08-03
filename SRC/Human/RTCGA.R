


source("https://bioconductor.org/biocLite.R")
#biocLite("RTCGA")
library(RTCGA)
infoTCGA()

biocLite(COAD.clinical)

#biocLite(RTCGA.mutations)
#biocLite(RTCGA.rnaseq)
#biocLite(RTCGA.CNV)
#biocLite(RTCGA.RPPA)
#biocLite(RTCGA.mRNA)
#biocLite(RTCGA.miRNASeq)
#biocLite(RTCGA.methylation)



library(RTCGA.clinical)
?clinical

clin <- survivalTCGA(COAD.clinical, READ.clinical,
                     extract.cols="admin.disease_code")
# Show the first few lines
head(clin)
table(clin$admin.disease_code)
xtabs(~admin.disease_code+patient.vital_status, data=clin) %>% addmargins()

coxph(Surv(times, patient.vital_status)~admin.disease_code, data=clin)
sfit <- survfit(Surv(times, patient.vital_status)~admin.disease_code, data=clin)
summary(sfit, times=seq(0,365*5,365))
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE)
#***************

```
```{r mutated genes}
bd.mutated <- read_delim("../data/frequently-mutated-genes.tsv", 
                         "\t", escape_double = FALSE, col_types = cols(`# Affected Cases Across the GDC` = col_number(), 
                                                                       `# Affected Cases in TCGA-COAD` = col_number()), 
                         trim_ws = TRUE)
bd.mutated$`# Affected Cases in TCGA-COAD` = bd.mutated$`# Affected Cases in TCGA-COAD` /400
bd.mutated$`# Affected Cases Across the GDC` = bd.mutated$`# Affected Cases Across the GDC` /10202

barplot(bd.mutated$`# Affected Cases in TCGA-COAD` , main=("Most Frequently mutated genes"),names.arg = bd.mutated$Symbol , col = 20 ,las = 2,xlim=c(0,40))



