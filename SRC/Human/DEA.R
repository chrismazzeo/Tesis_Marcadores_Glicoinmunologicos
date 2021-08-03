source("./src/functions/unzipHTCounts.r")
source("./src/functions/DEA_F.r")

#load raw counts matrix
rawCountMatrix = loadTSV(output.countMatrix.FileName)
rownames(rawCountMatrix) = rawCountMatrix$Ensembl_Gene_ID
rawCountMatrix = rawCountMatrix[,-1]

#load samples
samples = loadTSV(output.htseqCountsSamples)
samplesTypes = cbind(SampleID = samples$`Sample ID`, SampleType = samples$`Sample Type`)

#Run DE
deResults = DE_DESeq2 (rawCountMatrix, samplesTypes, ~SampleType, prefiltering = TRUE, coresNumber = 4, prefilteringCuttOff = settings.de.prefilteringCuttOff)
res = extractResults (dds, coresNumber = 4, pAdjustMethod = "BH", alpha = settings.de.alphaValue, contrast = c("SampleType","Primary Tumor","Solid Tissue Normal"),  format ="DataFrame",lfcThreshold = settings.de.lfcThreshold, filePath = output.DE.FileName)
summaryResults (res, filePath = settings.de.summaryFileName,padj=settings.de.alphaValue)

#filtering significant genes
res.filtered = significantpadjFiltering(res,padj= settings.de.alphaValue, filePath = output.DE.FilteredBySignificantFileName)

#Glyco filtering Homo Sampiens
HS_GlycoGenes = loadTSV(input.geneHSFilteringFileName)
res.glyco.HS.filtered = res[rownames(res) %in% HS_GlycoGenes$ensembl_gene_id,]
dim(res.glyco.HS.filtered)
saveTSV(as.data.frame(res.glyco.HS.filtered),output.DE.Glyco.HS.FileName)

#Glyco Filtering interseccion HS & MM
HS_MM_GlycoGenes = loadTSV(input.geneHS_MMFilteringFileName)
res.glyco.HS_MM.filtered = res[rownames(res) %in% HS_MM_GlycoGenes$HS_GENE_ID,]
dim(res.glyco.HS_MM.filtered)
saveTSV(as.data.frame(res.glyco.HS_MM.filtered),output.DE.Glyco.HS_MM.FileName)

#me falta guardar de los anterios los significativos solo para tenerlos 
#me falta guardar en la count matrix los genes que dieron significativos --> heatmap
#me falta guadar en la count matrix los genes de glyco que dieron significativos--> hetmap 
#me falta los heatmap by group