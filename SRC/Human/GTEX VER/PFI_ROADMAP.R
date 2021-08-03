library(readr)
library(Ecfun)
library(WGCNA)

#clear console
cat("\014") 
#clear enviroment
rm(list=ls())
list.of.packages <- c("readr","Ecfun","WGCNA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)




#source("http://bioconductor.org/biocLite.R") 
#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
#install.packages("WGCNA")



bd=  read_delim("./resources/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2", 
                "\t", escape_double = FALSE, trim_ws = TRUE)
bd_transposed <- transposeBigData(bd)
write_csv(bd_transposed,"./resources/Transposed_TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.csv")
#bd_transposed <- read.transpose("./resources/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2", header = TRUE,sep =  "\t")
View(bd_transposed[1:10,1:10])
dim(bd_transposed)




phenotype <- read_delim("./resources/TcgaTargetGTEX_phenotype.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
View(phenotype)
dim(phenotype)
print (unique(phenotype$detailed_category ))
print (unique(phenotype$`primary disease or tissue`))
print (unique(phenotype$`_sample_type`))




maxCluster = 47
clusters<-kmeans(bd_transposed, centers=maxCluster,iter.max = 100 )

#armar un archivo de settings y markdown

# tama?o de los Clusters
cat("\n********  Cluster size ******\n")
print(cbind(paste0("Cluster ",c(1:maxCluster),"  Size = ",clusters$size)))
#medias de los Clusters
print (cbind("centroides: ",c(1:maxCluster), round(clusters$centers,2)))
#grafico
#dev.new()

plot(data, col = clusters$cluster, pch = 19, main = paste0("K-means con k =",maxCluster))
points(clusters$centers, col = 1:maxCluster, pch = 8, cex = 5)

