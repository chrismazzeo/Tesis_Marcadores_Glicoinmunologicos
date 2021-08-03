#con este hacemos todos los heatmas lindos

library(readr)
library(writexl)
library(RColorBrewer)
library(pheatmap)

#dataset <- read_delim("/Volumes/Externo/Google Drive/Bioinformática/PFI/Tesina/Tesis Final/Resultados Finales/Consolidación/Tablas/FinalMergeLimpia_comparativa_ratones2.csv",  ";", escape_double = FALSE, trim_ws = TRUE)
#dataset <- read_delim("Results/Tablas/FinalMergeLimpia_humanos.csv", ";", escape_double = FALSE, trim_ws = TRUE)
dataset <- read_delim("~/Desktop/shares CCRAC/Todos los modelos-Tabla 1.csv", ";", escape_double = FALSE, trim_ws = TRUE)

View(dataset)


pdfOutput = "ratones_All_modelos.pdf"
dataFrameOutput= "resultRatonesModelos.csv"
minLog2FC = 2

                                    
dataset = as.data.frame(dataset)

heatmap.breakUp = 10
heatmap.breakDown = 7
heatmap.colorUp = colorRampPalette((brewer.pal(n = 8, name = "Reds")))(heatmap.breakUp)
heatmap.colorDown = colorRampPalette(rev(brewer.pal(n = 8, name = "Blues")))(heatmap.breakDown)
heatmap.colorZero = "white"
heatmap.colorNa = "white"
heatmap.colors=c( heatmap.colorDown,heatmap.colorZero,heatmap.colorUp)

a = dataset
if (filter){ #lo usamos para filtrar filas que no cumplan con el log2fc
    #removemos filas que no tengan lo2fC > minLog2FC
    a = dataset[,3:dim(dataset)[2]]
    
    
    
     heatmap.breaks = unique(c(
       seq(min(a, na.rm = TRUE),-1, length=heatmap.breakDown),
       0,
       rev(seq(max(a, na.rm = TRUE),1, length=heatmap.breakUp))
     ))
     a[abs(a) < minLog2FC] = NA
     a = dataset[rowSums(is.na(a)) != ncol(a), ] #remuevo filas con todo NA
     a[a==0] = NA
}
 
 
family.colors = settings.graphics.colors[3:(length(levels(as.factor(a$Family)))+2)]
names(family.colors) = levels(as.factor(a$Family))
annotationColor = list(
                       Family = family.colors)

annotationRow =data.frame(Family = a$Family)
rownames(annotationRow) = a$Genes
matrix = a[3:dim(a)[2]]
rownames(matrix) = a$Genes
#separación por familia de enzimas
rowGaps.freq = as.data.frame(table(a$Family))
rowGaps.freq = rowGaps.freq$Freq
rowGaps = 0
for (i in 1:length(rowGaps.freq)){
  rowGaps = rowGaps+ rowGaps.freq[i]
  rowGaps.freq[i] = rowGaps
}
rowGaps.freq = rowGaps.freq[1:(length(rowGaps.freq)-1)]
pheatmap(matrix, 
         cellheight = 10,cellwidth = 30,
         treeheight_row = 45, treeheight_col = 15, 
        annotation_row = annotationRow,
        annotation_colors = annotationColor,
        na_col ="white",
        display_numbers = FALSE, 
       # gaps_row = rowGaps.freq,
        fontsize_row = 8,
        fontsize_number = 4,
        cluster_rows = FALSE, cluster_cols = FALSE,  show_rownames = TRUE, show_colnames = TRUE,
        filename = paste0("./../",pdfOutput)
        )
HELPER_SAVE_DATA_FRAME(a,  paste0("./../",dataFrameOutput))

