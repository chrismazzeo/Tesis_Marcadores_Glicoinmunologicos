
mm = HELPER_LoadCSV( "./resources/GlycoGeneList/Glyco_MM_GeneList.csv")
mm$external_gene_name = toupper(mm$external_gene_name)
dim(mm)



geneAnnotationAttributes.immune = c('ensembl_gene_id','external_gene_name')
geneAnnotation.immune = ANN_GetAnnotationAll(specie = prject.rnaseq.ensemblSpecie, attributes = geneAnnotationAttributes.immune,collapseInRow = FALSE)

new = geneAnnotation.immune[geneAnnotation.immune$external_gene_name %in% toupper(mm$external_gene_name),]

a = merge(new, mm,by.x = 'external_gene_name', by.y = 'external_gene_name')

a = data.frame (ensembl_gene_id = a$ensembl_gene_id.x, external_gene_name = a$external_gene_name, Group = a$Group )
HELPER_SaveCSV(a, "./caca.csv", rowNames = FALSE)



fpkm = HELPER_LoadCSV(samplesListFPKMPath)




library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.mRNA)
library(RTCGA.rnaseq)




#ver esto
#sto para una maezcla de go y gsea
library(GOplot) #https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html
data(EC)
head(EC$david)
head(EC$genelist)
circ <- circle_dat(EC$david, EC$genelist)
GOBar(subset(circ, category == 'BP'))
GOBar(circ, display = 'multiple')
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))
GOBubble(circ, labels = 3)
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)




VERRRR

#**********************
#Reactome
#**********************
#TODO falta terminar
reactomeList = c('ensembl_gene_id', 'entrezgene')
reactomeList.result = ANN_GetAnnotation(specie = "mmusculus_gene_ensembl",filterIDs = myInterestedGenes, ensemblIDType = "ensembl_gene_id",attributes = reactomeList, uniqueRows = TRUE)
reactomeList.result = reactomeList.result[!duplicated(reactomeList.result$entrezgene),]
reactomeList.result.stat = merge(reactomeList.result, res.annotation, by = "ensembl_gene_id")
reactomeList.result.stat = data.frame(reactomeList.result.stat$entrezgene, reactomeList.result.stat$stat)

#ReactomePA
library(ReactomePA)
yy = enrichPathway(reactomeList.result.stat$reactomeList.result.stat.entrezgene, organism = "mouse", pvalueCutoff = 0.05,readable = TRUE)
View(summary(yy))

reactomeList.result.stat.ordered = reactomeList.result.stat[order(reactomeList.result.stat$reactomeList.result.stat.stat, decreasing = TRUE),]

xx = gsePathway(reactomeList.result.stat.ordered$reactomeList.result.stat.entrezgene, organism = "mouse")

viewPathway(, organism = "mouse")









#************************************






onts = c( "MF", "BP", "CC" )

## prepare data
sampleGOdata <- new( "topGOdata", ontology="CC", allGenes = geneList, nodeSize=10,
                     annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl" )

#clasig test
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "Fisher" )

#Kolmogorov-Smirnov test
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")




## look at results
allRes <- GenTable( sampleGOdata,
                    classicFisher = resultFisher,
                    classicKS = resultKS, 
                    elimKS = resultKS.elim,
                    
                    orderBy = "elimKS" , ranksOf = "classicFisher", topNodes = 1500)
View(allRes)

pValue.classic <- score(resultKS)


pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
      elim = pValue.elim[sel.go],
      classic = pValue.classic[sel.go])






