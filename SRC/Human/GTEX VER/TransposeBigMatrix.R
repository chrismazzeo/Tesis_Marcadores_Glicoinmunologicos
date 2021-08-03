#clear console
cat("\014") 
#clear enviroment
rm(list=ls())

list.of.packages <- c("readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(readr)


fileName = "TCGA_COAD_HTCounts.csv"
sourcePath = "Results/"

filePath = paste0(sourcePath,fileName)


bd=  read_delim(filePath,",", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)

bd_transposed <- t(bd)
columns = bd_transposed[1,]
bd_transposed = bd_transposed[-1,]
colnames(bd_transposed) = columns



write.csv(bd_transposed, file =  paste0(sourcePath,"Transposed_",fileName))


dim(bd_transposed)




