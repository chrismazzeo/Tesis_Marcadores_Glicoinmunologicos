
bioconductor.packages <- c("topGO","goseq","DEXSeq","tximport","pathview", "gage", "gageData", "GenomicAlignments","rhdf5","limma","GO.db",
                           "fgsea","GOplot","RDAVIDWebService","ReactomePA","PGSEA","gskb","IHW","apeglm","clusterProfiler",
                           "Biobase","affy","gcrma","limma","GEOquery","affyPLM","ArrayExpress","genefilter","geneplotter","ReactomePA",
                           "GSVA",
                           "TxDb.Hsapiens.UCSC.hg19.knownGene","TxDb.Mmusculus.UCSC.mm10.knownGene",
                           "EnhancedVolcano",
                           "AnnotationDbi",
                           "DBI",
                           "frma",
                           "org.Mm.eg.db",
                           "org.Hs.eg.db",
                           "mouse4302.db","moe430a.db",
                           "mouse4302mmentrezg.db","mouse4302mmentrezgcdf", "mouse4302mmentrezgprobe","mouse4302frmavecs",
                           "DaMiRseq",
                           "sSeq",
                           "Pigengene","biosigner","DEXUS","pcaGoPromoter"


                           )
# ,"RTCGA","RTCGA.clinical","RTCGA.mRNA","RTCGA.rnaseq", "RTCGA.sample"

new.packages <- bioconductor.packages[!(bioconductor.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  BiocManager::install(new.packages, version = "3.9",ask = FALSE)
}

#sudo ln -f -s $(/usr/libexec/java_home)/lib/server/libjvm.dylib /usr/local/lib  para rJava
list.of.packages <- c("RColorBrewer","pheatmap","pca3d","factoextra","FactoMineR","rgl","Rtsne","plotly","NbClust","BiocManager","dplyr","tidyverse","VennDiagram","pkgdown","grid",
                      "psych","corrplot","asbio","mvShapiroTest","ggpubr","kableExtra","FSA","gridExtra","ggsci","fpc","cluster","readr","reshape2","mvnormtest","Rserve","e1071","parallel","preprocessCore","colorRamps",
                      "survival","survminer", "rgl","writexl","NBLDA","tree")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)




#TODO agregar los que faltan
library(BiocManager)
library(RColorBrewer)
library(pheatmap)
library(pca3d)
library(factoextra)
library(FactoMineR)
library(rgl)
library(Rtsne)
library(plotly)
library(NbClust)
library(tximport)
library(DEXSeq)
library(topGO)
library(goseq)
library(GOplot)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GO.db)
library(dplyr)
library(tidyverse)
library(apeglm)
library(limma)
library(VennDiagram)
library(fgsea)
library(clusterProfiler)
library(pathview)
library(AnnotationDbi)
library(DBI)
library(readr)
library(cluster)
library(fpc)
library(ggsci)
library(grid)
library(gridExtra)


library (FSA)
library(kableExtra)
library(ggpubr)
library(mvShapiroTest)
library(asbio)
library(corrplot)
library(psych)
library(reshape2)
library(mvnormtest)

library(GEOquery)
library(limma)
library(Biobase)
library(affy)
library(gcrma)
library(affyPLM)
library(moe430a.db)
library(mouse4302.db)
library(frma)

library("IHW")
library(ArrayExpress)
library(genefilter)
library(geneplotter)
#library(ReactomePA)
library(GSVA)
library(pkgdown)
library(survival)
library(survminer)
library(rgl)
library(car)
library(writexl)
library(Rfast)
library(EnhancedVolcano)


source("./src/shared/functions/DEA_Function.r")
source("./src/shared/functions/DA_Function.r")
source("./src/shared/functions/HELPER_Function.r")
source("./src/shared/functions/ANNOTATIONS_Function.r")
source("./src/shared/functions/FILTERING_Function.r")
source("./src/shared/functions/IMMUNE_Function.r")
source("./src/shared/functions/ONTHO_Function.r")
source("./src/shared/functions/SURVIVAL_FUNC.r")
source("./src/shared/functions/GRAPHICS_Function.r")
source("./src/shared/functions/FPKM_TPM_Function.r")
source("./src/shared/functions/DEA_CONSOLIDADO_Function.r")
source("./ExternalLibs/ImmuCC/microarray/Microarray_Deconvolution.r")
source("./ExternalLibs/CIBERSORT/CIBERSORT.r")



