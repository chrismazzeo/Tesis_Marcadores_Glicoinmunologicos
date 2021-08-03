#***************************************************************
#Clinical
#***************************************************************


gdc_sample <- read_delim("~/Desktop/gdc_sample_sheet.2018-08-17.tsv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
clinical <- read_delim("~/Desktop/clinical.cart.2018-08-17/clinical.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

exposure <- read_delim("~/Desktop/clinical.cart.2018-08-17/exposure.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

#***************************************************************
#Biospecimen
#***************************************************************
slide <- read_delim("/Volumes/Externo/Projects/R pipelines/PFI ROADMAP/Resources/GDC/TCGA/COAD/biospecimen/slide.tsv", 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
# samples <- read_delim("/Volumes/Externo/Projects/R pipelines/PFI ROADMAP/Resources/GDC/TCGA/COAD/biospecimen/sample.tsv", 
#                      "\t", escape_double = FALSE, trim_ws = TRUE)
# 
# portion <- read_delim("/Volumes/Externo/Projects/R pipelines/PFI ROADMAP/Resources/GDC/TCGA/COAD/biospecimen/portion.tsv", 
#                       "\t", escape_double = FALSE, trim_ws = TRUE)
# analyte <- read_delim("/Volumes/Externo/Projects/R pipelines/PFI ROADMAP/Resources/GDC/TCGA/COAD/biospecimen/analyte.tsv", 
#                       "\t", escape_double = FALSE, trim_ws = TRUE)
# 
# aliquot <- read_delim("/Volumes/Externo/Projects/R pipelines/PFI ROADMAP/Resources/GDC/TCGA/COAD/biospecimen/aliquot.tsv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)

#aca hay repeplicas que hay que volar para facilitar todo
samples = read_delim(input.htSeqSampleFileName,"\t", escape_double = FALSE, trim_ws = TRUE)


infiltration = as.data.frame(slide[slide$sample_submitter_id %in% samples$`Sample ID`,])
infiltration = cbind(infiltration, Group =c(1:dim(infiltration)[1]) )
#add group
for (i in 1:dim(infiltration)[1]){
    infiltration$Group[i] = samples$`Sample Type`[samples$`Sample ID` == infiltration$sample_submitter_id[i]] 
}
#
# infiltration = slides filtered + group
#
a esto le tengo que pasar las cosas y que me devuelva la matrix,  de ahi llamo a otras funciones para graficar


que hago con esto?
  
  
biospecimenInfiltration = function (infiltration){
  tumorPurity = data.frame(Normal = as.numeric(infiltration$percent_normal_cells), Stromal = as.numeric(infiltration$percent_stromal_cells) ,Tumor = as.numeric(infiltration$percent_tumor_cells),Necrosis = as.numeric(infiltration$percent_necrosis),  InflamatoryInfiltration = as.numeric(infiltration$percent_inflam_infiltration))
  tumorPurity = cbind(SampleID =infiltration$sample_submitter_id, Group = infiltration$Group, SectionLocation = infiltration$section_location, tumorPurity, "Total (%)" =rowSums(tumorPurity, na.rm = TRUE))
  tumorPurity = tumorPurity[order(tumorPurity$Group, tumorPurity$SampleID,tumorPurity$SectionLocation, decreasing = FALSE),]
  barplot (height=as.matrix(t(tumorPurity[,4:7])), col = settings.color[1:4])
  
  infiltrationMatrix = data.frame (Monocyte = as.numeric(infiltration$percent_monocyte_infiltration),Eosinophil = as.numeric(infiltration$percent_eosinophil_infiltration), Lymphocyte = as.numeric(infiltration$percent_lymphocyte_infiltration), Neutrophil = as.numeric(infiltration$percent_neutrophil_infiltration),Granulocyte =  as.numeric(infiltration$percent_granulocyte_infiltration))
  infiltrationMatrix = cbind(SampleID =infiltration$sample_submitter_id, Group = infiltration$Group,SectionLocation = infiltration$section_location,infiltrationMatrix, "Total (%)" =rowSums(infiltrationMatrix, na.rm = TRUE))
  infiltrationMatrix = infiltrationMatrix[order(infiltrationMatrix$Group, infiltrationMatrix$SampleID,infiltrationMatrix$SectionLocation, decreasing = FALSE),]
}








#sample Type
group = sample$sample_type[(sample$sample_type == "Primary Tumor" | sample$sample_type == "Solid Tissue Normal")]
tumorPurity.filtered = cbind (tumorPurity, SampleType = group)
infiltrationMatrix.filtered = cbind (infiltrationMatrix, SampleType = group)

#filtrado x sample type
tumorPurity.filtered = tumorPurity[(sample$sample_type == "Primary Tumor" | sample$sample_type == "Solid Tissue Normal"),]
infiltrationMatrix.filtered = infiltrationMatrix[(sample$sample_type == "Primary Tumor" | sample$sample_type == "Solid Tissue Normal"),]

#filtrado total >0
tumorPurity.filtered  = tumorPurity.filtered[tumorPurity.filtered$Total>0,]
infiltrationMatrix.filtered  = infiltrationMatrix.filtered[infiltrationMatrix.filtered$Total>0,]

View(tumorPurity.filtered)
View(infiltrationMatrix.filtered)
barplot (inflamatoryCellMatrix)


