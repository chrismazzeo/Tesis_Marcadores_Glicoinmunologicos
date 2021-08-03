CreateGroupFolder = function (outputDir, groupTest){
  
  dir.create(file.path(outputDir), showWarnings = FALSE)
  dir.create(file.path(paste0(outputDir,"/CountMatrix/")), showWarnings = FALSE)
  dir.create(file.path(paste0(outputDir,"/DescriptiveAnalysis/")), showWarnings = FALSE)
  dir.create(file.path(paste0(outputDir,"/DEA/")), showWarnings = FALSE)
  
  for (i in 1:length(groupTest)){
    group = paste0(groupTest[[i]][2],"_VS_",groupTest[[i]][3])
    outputDirGroup = paste0(outputDir,"/DEA/",group,"/" )
    dir.create(file.path(outputDirGroup), showWarnings = FALSE)
    
    dir.create(file.path(paste0(outputDirGroup,"DEA/")), showWarnings = FALSE)
    dir.create(file.path(paste0(outputDirGroup,"DEA/GlycoHeatmaps/")), showWarnings = FALSE)
    dir.create(file.path(paste0(outputDirGroup,"DEA/GlycoHeatmaps/Log2FC_1/")), showWarnings = FALSE)
    
    dir.create(file.path(paste0(outputDirGroup,"GO/")), showWarnings = FALSE)
    dir.create(file.path( paste0(outputDirGroup,"GSEA/")), showWarnings = FALSE)
    dir.create(file.path(paste0(outputDirGroup,"KEGG/")), showWarnings = FALSE)
    
    if  (project.settings.DEA.H0.Log2FC >0){
        outputDirGroup = paste0(outputDir,"/DEA/",group,"_Log2FC_",project.settings.DEA.H0.Log2FC,"/")
        dir.create(file.path(outputDirGroup), showWarnings = FALSE)
        
        dir.create(file.path(paste0(outputDirGroup,"DEA/")), showWarnings = FALSE)
        dir.create(file.path(paste0(outputDirGroup,"DEA/GlycoHeatmaps/")), showWarnings = FALSE)
        dir.create(file.path(paste0(outputDirGroup,"DEA/GlycoHeatmaps/Log2FC_1/")), showWarnings = FALSE)
        
        dir.create(file.path(paste0(outputDirGroup,"GO/")), showWarnings = FALSE)
        dir.create(file.path( paste0(outputDirGroup,"GSEA/")), showWarnings = FALSE)
        dir.create(file.path(paste0(outputDirGroup,"KEGG/")), showWarnings = FALSE)
    }
  }
}

CreateInitialFolder = function (specie,projectName,resourcePath, resultPath ){
  dir.create(file.path( resultPath), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path( paste0(resultPath,projectName)), showWarnings = FALSE)
  
  dir.create(file.path( paste0(resultPath,projectName,"/Samples/")), showWarnings = FALSE)
  
  return (paste0(resultPath,projectName,"/"))
}

project.resultsPath = CreateInitialFolder(specie = project.specie,
                                          projectName = project.name,
                                          resultPath = project.resultsPath)

#*******************************************************************************************************
# Samples
#*******************************************************************************************************
geneFilteringList =  HELPER_LoadCSV(geneListFilterPath)
HELPER_SAVE_DATA_FRAME(geneFilteringList, paste0(project.resultsPath,"GeneFilteringList.csv")) 

#*******************************************************************************************************
# Gene Filtering List
#*******************************************************************************************************
HELPER_SAVE_DATA_FRAME(samplesList, paste0(project.resultsPath,"/Samples/",project.name,"_SampleList.csv"), rowNames = FALSE) 

#*******************************************************************************************************
# Annotations
#*******************************************************************************************************
#geneLevel
if (project.rnaseq.runHTSEQ | prject.rnaseq.runKallisto | project.rnaseq.runDEXSeq) {
    geneAnnotationAttributes = c('ensembl_gene_id_version','ensembl_gene_id','entrezgene','external_gene_name','gene_biotype','description','chromosome_name','start_position','end_position','strand')
    geneAnnotation = ANN_GetAnnotationAll(specie = prject.rnaseq.ensemblSpecie, attributes = geneAnnotationAttributes,collapseInRow = TRUE)
    
    #transcriptLevel
    transcriptAnnotationAttributes = c('ensembl_transcript_id_version','ensembl_gene_id')
    transcriptAnnotation = ANN_GetAnnotationAll(specie = prject.rnaseq.ensemblSpecie, attributes = transcriptAnnotationAttributes,collapseInRow = FALSE)
    
    if (project.specie != "HUMAN"){
      #mouse-human homologs
      homologAnnotationAttributes = c('ensembl_gene_id','hsapiens_homolog_associated_gene_name')
      homologAnnotation = ANN_GetAnnotationAll(specie = prject.rnaseq.ensemblSpecie, attributes = homologAnnotationAttributes,  collapseInRow = FALSE)
    }
}

if (project.type == "RNASEQ"){
    if (project.rnaseq.runHTSEQ){
      source("./src/shared/rnaseq/HTSEQ.r")
    }
     if (prject.rnaseq.runKallisto){
       source("./src/shared/rnaseq/KALLISTO.r")
     }
    if (project.rnaseq.runDEXSeq){
      source("./src/shared/rnaseq/DEXSEQ.r")
    }
    if (project.rnaseq.run.immune){
        if  (project.specie == "MOUSE"){
          source("./src/shared/rnaseq/MM_IMMUCC.r")
        }
        if (project.specie == "HUMAN"){
          source("./src/shared/rnaseq/HS_RNA_SEQ_ImmuneInfiltrations.r")
        }
    }
}
if (project.type == "TCGA"){
  if (project.rnaseq.runHTSEQ){
     source("./src/shared/tcga/tcga.r")
  }
   if(project.rnaseq.run.immune){
      source("./src/shared/rnaseq/HS_RNA_SEQ_ImmuneInfiltrations.r")
   }

}