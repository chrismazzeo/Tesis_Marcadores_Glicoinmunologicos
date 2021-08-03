install.packages("PerformanceAnalytics")
install.packages("corrplot")
install.packages("Hmisc")
library("PerformanceAnalytics")
library(corrplot)
library("Hmisc")
library(rstatix)
library(psych)
library(readr)

db <- read_delim("~/Desktop/Raton-HumanoFinalMergeTesis.csv", 
                                          ";", escape_double = FALSE, trim_ws = TRUE)


db = db[,-1:-2]
db = db[,-11]
View(db)

#Test normality
mshapiro_test(db)
#Calculate correlation with pvalue
db.res <- rcorr(as.matrix(db), type ="spearman")
db.table = flattenCorrMatrix(db.res$r, db.res$P)
#save table
View(db.table)
HELPER_SAVE_DATA_FRAME(as.data.frame(db.table), "~/Desktop/Correlation.csv")

#correlation plot1
pdf("~/Desktop/Correlation_chart.pdf", width = 20, height= 10)
chart.Correlation(db, histogram=TRUE, pch=1,method = "spearman")
dev.off()

#correlation plot 2
pdf("~/Desktop/Correlation_chart.pdf", width = 20, height= 10)
pairs.panels(db, scale=TRUE, method="spearman")
dev.off()
#correlation plot 3
settings.heatmap.correlation.color=colorRampPalette(c("blue","white","brown"))(10)

pdf("~/Desktop/Correlation.pdf", width = 10, height= 12)
corrplot(corr = db.res$r,
          p.mat = db.res$p,
          type = "lower",
          sig.level = 0.05,
          insig = "blank",
          title ="Correlation", 
         diag = FALSE,
          col = settings.heatmap.correlation.color,
          tl.col="black", #Text label color and rotation
          number.cex =0.4,
          cl.cex = 1,
          na.label = "  ",
          mar=c(0,0,1,0)
 ) 
 dev.off()
 



# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
 ut <- upper.tri(cormat)
 data.frame(
   row = rownames(cormat)[row(cormat)[ut]],
   column = rownames(cormat)[col(cormat)[ut]],
   cor  =(cormat)[ut],
   p = pmat[ut]
 )
}
