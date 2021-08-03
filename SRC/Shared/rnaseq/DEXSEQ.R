#*******************************************************************************************************
#TODO DExsEQ independiente de la especie

#*******************************************************************************************************
#splicing alternativo
#*******************************************************************************************************
dir.create(file.path(paste0(project.resultsPath,"DEXSeq")), showWarnings = FALSE)
dir.create(file.path(paste0(project.resultsPath,"DEXSeq/plots/")), showWarnings = FALSE)
dir.create(file.path(paste0(project.resultsPath,"DEXSeq/DEA")), showWarnings = FALSE)

dexseq.res =  DEXSEQ_CALCULATE(
  files =  paste0(project.resourcePath,"dexseq/",samplesList$HTSEQ_FILE),
  sampleList = samplesList,
  flattenedfile = "./resources/mouse/Mus_musculus.GRCm38.94.gtf.DEXSeq.chr.gff",
  design =  ~ SAMPLE  + exon + GROUP:exon,
  fitExpToVar =  "GROUP",
  coresNumber = setttings.coresNumber
)


dexeq.significantID = HELPER_FILER_SAVE_DEXEQ_PFI(dexseq.res,geneFilteringList$ensembl_gene_id, padj = project.settings.DEA.padj, outputDir = paste0(project.resultsPath,"DEXSeq/"))
DEXSEQ_PLOT_ALL_PDF(dexseq.res,dexeq.significantID$significantID.all,fitExpToVar= "GROUP",padj = project.settings.DEA.padj, outputDir = paste0(project.resultsPath,"DEXSeq/plots/")) 


