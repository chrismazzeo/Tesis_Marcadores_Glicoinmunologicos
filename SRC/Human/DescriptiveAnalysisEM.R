source("./src/functions/unzipHTCounts.r")
source("./src/functions/removerOutliersFunc.r")
#install.packages("NbClust")
library(pheatmap)
library(pca3d)
library(factoextra)
library(rgl)
library(Rtsne)
library(plotly)
library(NbClust)
#******************************************************************************************************************************
#Load samples & count matrix
#******************************************************************************************************************************
samples = loadTSV(output.htseqCountsSamples)
dim(samples)
countMatrix = loadTSV(output.countMatrix.FileName)
genesID = countMatrix[,1]
#removemos los Genes IDs
countMatrix = countMatrix[,-1]

#transpose
countMatrix =t(countMatrix) #row = samples
#******************************************************************************************************************************
#removemos samples outliers con cluster de 1 o 2  , esto a ojo
outliersList = c("TCGA-A6-5656-01A","TCGA-A6-3809-01A", "TCGA-A6-3810-01A")
outliersList2 = c("TCGA-A6-5656","TCGA-A6-3809", "TCGA-A6-3810")
countMatrix = countMatrix[!(rownames(countMatrix) %in% outliersList),]
samples = samples[!(samples$`Case ID` %in% outliersList2),]
dim(countMatrix)
dim(samples)

#******************************************************************************************************************************
#filter colums = 0
countMatrix = countMatrix[,colSums(countMatrix)>settings.de.prefilteringCuttOff]
dim(countMatrix)
#******************************************************************************************************************************
#scale
countMatrix = scale(countMatrix)
#******************************************************************************************************************************
#tSNE
#******************************************************************************************************************************
set.seed(123)
#2d
tsne.2d <- Rtsne(countMatrix, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)

optimalCluster.elbow = fviz_nbclust(tsne.2d$Y, kmeans, method = "wss") +
  labs(subtitle = "Elbow method")
optimalCluster.silhoutte = fviz_nbclust(tsne.2d$Y, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")
optimalCluster.gap_stat = fviz_nbclust(tsne.2d$Y, kmeans, method = "gap_stat",nboot = 50) +
  labs(subtitle = "Gap statistic method")
nbClust <- NbClust(tsne.2d$Y, distance = "euclidean", min.nc = 2,max.nc = 10, method = "kmeans")
optimalCluster.nbClust = fviz_nbclust(nbClust)+ labs(subtitle = "NbClust")

pdf(output.tsne.2d.optimalClusterFileName)
optimalCluster.elbow
optimalCluster.silhoutte
optimalCluster.gap_stat
optimalCluster.nbClust
dev.off()

cluster.tsne2d<-kmeans(scale(tsne.2d$Y), centers=3,iter.max = settings.kmeans.maxIteration, nstart = 20)
table(cluster.tsne2d$cluster)
tsne.graphics.2d = plot_ly(as.data.frame(tsne.2d$Y), x = ~V1, y = ~V2, color = settings.color[cluster.tsne2d$cluster])
orca(tsne.graphics.2d, output.tsne.2d.resultFileName, format = "pdf")

#3d
tsne.3d <- Rtsne(countMatrix, dims = 3, perplexity=30, verbose=TRUE, max_iter = 500)

optimalCluster.elbow = fviz_nbclust(tsne.3d$Y, kmeans, method = "wss")+
  labs(subtitle = "Elbow method")
optimalCluster.silhoutte = fviz_nbclust(tsne.3d$Y, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")
optimalCluster.gap_stat = fviz_nbclust(tsne.3d$Y, kmeans, method = "gap_stat",nboot = 50) +
  labs(subtitle = "Gap statistic method")
nbClust <- NbClust(tsne.3d$Y, distance = "euclidean", min.nc = 2,max.nc = 10, method = "kmeans")
optimalCluster.nbClust = fviz_nbclust(nbClust)+ labs(subtitle = "NbClust")

pdf(output.tsne.3d.optimalClusterFileName)
optimalCluster.elbow
optimalCluster.silhoutte
optimalCluster.gap_stat
optimalCluster.nbClust
dev.off()

cluster.tsne3d<-kmeans(scale(tsne.3d$Y), centers=3,iter.max = settings.kmeans.maxIteration, nstart = 20)
table(cluster.tsne3d$cluster)
tsne.graphics.3d = plot_ly(as.data.frame(tsne.3d$Y), x = ~V1, y = ~V2, z = ~V3, color = settings.color[cluster.tsne3d$cluster])
orca(tsne.graphics.3d, output.tsne.3d.resultFileName, format = "pdf")


#******************************************************************************************************************************
#sample-to-sample distance  muestra que samples estan mas cercanas
#******************************************************************************************************************************
sampleDists <-  as.matrix(dist(countMatrix))
rownames(sampleDists) <- samples$`Sample Type`
colnames(sampleDists) <- NULL
pheatmap(sampleDists,  
         main = "Sample Distance\n",
         fontsize_row =1,
         
         col=settings.heatmp.distance.colorRamp,
         filename = output.heatmap.sampleDistanceFileName)
#******************************************************************************************************************************
#Kmeans
#******************************************************************************************************************************
#cantidad de clusters Elbow method
optimalCluster.elbow = fviz_nbclust(countMatrix, kmeans, method = "wss")+
  labs(subtitle = "Elbow method")
#optimalCluster.silhoutte = fviz_nbclust(countMatrix, kmeans, method = "silhouette") +  labs(subtitle = "Silhouette method")
#optimalCluster.gap_stat = fviz_nbclust(countMatrix, kmeans, method = "gap_stat",nboot = 50) +   labs(subtitle = "Gap statistic method")
nbClust <- NbClust(countMatrix, distance = "euclidean", min.nc = 2,max.nc = 4, method = "kmeans")
optimalCluster.nbClust = fviz_nbclust(nbClust)+ labs(subtitle = "NbClust")

pdf(output.kmeans.optimalClustersFileName)
optimalCluster.elbow
#optimalCluster.silhoutte
#optimalCluster.gap_stat
#optimalCluster.nbClust
dev.off()

#kmeans
clusters<-kmeans(countMatrix, centers=3,iter.max = settings.kmeans.maxIteration, nstart = 20)
table(clusters$cluster)
objectsInClusters = rbind(paste0("Cluster ",c(1:settings.kmeans.maxCluster),"  Size = ",clusters$size))

#ploteamos cluster
pdf(output.kmeans.resultFileName)
plot(countMatrix, col = settings.color[1:length(clusters$size)], pch = 19, main = paste0("K-means  (k = ",settings.kmeans.maxCluster,")"))
legend("topright",objectsInClusters, cex = 0.8,fill = settings.color[1:length(clusters$size)])
dev.off()

#******************************************************************************************************************************
#PCA
#******************************************************************************************************************************
pca <- prcomp(countMatrix, center = TRUE)
pdf(output.pca.PCPorcentage)
fviz_eig(pca, addlabels=TRUE, hjust = -0.3, title = "Porcentage explained for each PC")+ theme_minimal()
dev.off()
type = samples$`Sample Type`
type[type == "Primary Tumor"] = 1
type[type == "Solid Tissue Normal"] = 2
type = as.numeric(type)
pca3d(pca,new = TRUE, title = "TCGA COAD 3D PCA", show.ellipses = TRUE, show.axe.titles = TRUE, show.plane = FALSE, group =samples$`Sample Type` , 
      col = settings.color[type],
      legend =legend3d("center", c("2D Points", "3D Points")))
snapshotPCA3d(output.pca.3dFileName)


pdf(output.pca.2dFileName)
pca2d(pca,new = TRUE, title = paste0("TCGA COAD 2D PCA"), group =samples$`Sample Type`,  
      col = settings.color[type],show.ellipses = FALSE,  show.axe.titles = TRUE, show.plane = FALSE)
legend("topright",levels(factor(samples$`Sample Type`)),   col = settings.color[1:length(type)], fill=  settings.color[1:length(type)])
dev.off()

#***************************************************************
#guardamos los clusters
#***************************************************************
cluster.result = data.frame(sample = rownames(countMatrix),type = samples$`Sample Type`, KmeanCluster =clusters$cluster, TSNECluster2d = cluster.tsne2d$cluster,TSNECluster3d = cluster.tsne3d$cluster )
saveTSV(cluster.result,output.kmeans.custerFileName)

#******************************************************************************************************************************
#PCA + Kmeans
#******************************************************************************************************************************
pca3d(pca,new = FALSE, title = "3D PCA" ,col = settings.color[clusters$cluster],show.ellipses = TRUE,show.group.labels = FALSE, show.axe.titles = TRUE)
snapshotPCA3d(paste0(folder.output.da,project.folder,"_PCA_3DWithKmeans.png"))
dev.off()



