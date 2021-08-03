library(tidyverse)
library(stringr)

source("./src/functions/unzipHTCounts.r")
source("./src/functions/FPKM_TPM.r")
source("./src/functions/AnnotationsExtra.r")
source("./src/functions/Infiltration.R")

#********************************************
#0. merge clinical data
#********************************************

clinical = loadTSV(input.clinical)
diagnoses = loadTSV(input.clinical.exposure)
diagnoses = as.data.frame(cbind(submitter_id = diagnoses$submitter_id,weight = diagnoses$weight,height = diagnoses$height,bmi = diagnoses$bmi))
clinical = merge (x = clinical,
                   y = diagnoses,
                   by.x = "submitter_id",
                   by.y ="submitter_id" )



All_CDEs <- loadTSV(input.clinical.msistatus)
All_CDEs = as.data.frame(All_CDEs)
rownames(All_CDEs) = All_CDEs[,1]
All_CDEs = All_CDEs[,-1]
All_CDEs =as.data.frame(t(All_CDEs))
rownames(All_CDEs) = toupper(rownames(All_CDEs))
All_CDEs = cbind (submitter_id = rownames(All_CDEs),All_CDEs)


msi =  data.frame(submitter_id = All_CDEs$submitter_id,All_CDEs$anatomic_neoplasm_subdivision, All_CDEs$circumferential_resection_margin,
                           All_CDEs$braf_gene_analysis_performed,All_CDEs$braf_gene_analysis_result,All_CDEs$colon_polyps_present,
                           All_CDEs$kras_gene_analysis_performed,All_CDEs$kras_mutation_codon,All_CDEs$kras_mutation_found,All_CDEs$loss_expression_of_mismatch_repair_proteins_by_ihc,
                           All_CDEs$loss_expression_of_mismatch_repair_proteins_by_ihc_result,All_CDEs$`loss_expression_of_mismatch_repair_proteins_by_ihc_result-2`,
                           All_CDEs$`loss_expression_of_mismatch_repair_proteins_by_ihc_result-3`,All_CDEs$`loss_expression_of_mismatch_repair_proteins_by_ihc_result-4`,
                           All_CDEs$lymph_node_examined_count,All_CDEs$lymphatic_invasion,All_CDEs$primary_lymph_node_presentation_assessment,All_CDEs$number_of_lymphnodes_positive_by_he,
                           All_CDEs$number_of_lymphnodes_positive_by_ihc,All_CDEs$microsatellite_instability,All_CDEs$mononucleotide_and_dinucleotide_marker_panel_analysis_status,
                           All_CDEs$mononucleotide_marker_panel_analysis_status,All_CDEs$non_nodal_tumor_deposits,All_CDEs$number_of_abnormal_loci,
                           All_CDEs$number_of_first_degree_relatives_with_cancer_diagnosis,All_CDEs$other_dx,All_CDEs$pathologic_m,
                           All_CDEs$pathologic_n,All_CDEs$pathologic_t,All_CDEs$pathologic_stage,All_CDEs$perineural_invasion_present,All_CDEs$person_neoplasm_cancer_status,
                           All_CDEs$postoperative_rx_tx,All_CDEs$preoperative_pretreatment_cea_level,All_CDEs$primary_therapy_outcome_success,All_CDEs$radiation_therapy,
                           All_CDEs$residual_tumor,All_CDEs$synchronous_colon_cancer_present,All_CDEs$tissue_source_site,All_CDEs$venous_invasion,All_CDEs$withdrawn,
                   stringsAsFactors=FALSE)
names = colnames(msi)
colnames(msi) = gsub("All_CDEs.","",names)
clinical = merge (x = clinical,
                   y = msi,
                   by.x = "submitter_id",
                   by.y ="submitter_id" ,  all.x = TRUE)

#**********************
#1. Samples Filtering
#**********************
htSeqCountsSamples = filterSamples(input.reads.htSeq.countsFile,output.htseqCountsSamples)
htSeqFPKMSamples = filterSamples(input.reads.htSeq.FPKMFile,output.htseqFPKMSamples)


clinical = clinical[clinical$submitter_id %in% htSeqCountsSamples$`Case ID`,]
saveTSV(clinical,output.clinicalSamples )

#********************************************
#2. Create Count Matrix for DE
#********************************************
htseqCountMatrix = createHTCountsMatrix(folder.input.reads.htSeq.countsFile, htSeqCountsSamples,output.countMatrix.FileName)  
htseqFPKMMatrix = createHTCountsMatrix(folder.input.reads.htSeq.FPKMFile, htSeqFPKMSamples,output.countMatrix.FPKM.FileName)

htseqTPMMatrix = htseqFPKMMatrix [,-1]
htseqTPMMatrix=  rpkmToTpm(htseqTPMMatrix)
htseqTPMMatrix = cbind (Ensembl_Gene_ID = htseqCountMatrix$Ensembl_Gene_ID,htseqTPMMatrix)
saveTSV(htseqTPMMatrix,output.countMatrix.TPM.FileName)
#********************************************
#3. Get  Annotations
#********************************************
annotationIDList = c('hgnc_symbol','gene_biotype','superfamily','family','family_description','phenotype_description','chromosome_name','start_position','end_position','description','strand','band')
annResult = getAnnotation(specie = "hsapiens_gene_ensembl",ensemblIDs= htseqCountMatrix$Ensembl_Gene_ID, ensemblIDType = "ensembl_gene_id",annotationIDList=annotationIDList, oneRow = TRUE)
annResult = collapseInRow(annResult)

annotationGeneOrthology = c('mmusculus_homolog_ensembl_gene','mmusculus_homolog_associated_gene_name','mmusculus_homolog_orthology_type')
annOrthologyResult = getAnnotation(specie = "hsapiens_gene_ensembl",ensemblIDs= htseqCountMatrix$Ensembl_Gene_ID, ensemblIDType = "ensembl_gene_id",annotationIDList=annotationGeneOrthology, oneRow = TRUE)
#annOrthologyResult = collapseInRow(annOrthologyResult)

saveTSV(annResult,output.annotations.hs)
saveTSV(annOrthologyResult,output.annotations.mouse)

#********************************************
#4. Create Matrix For Immune Infiltration
#********************************************

htseqCountMatrix.immune = createMatrixForImmuneInfiltration(htseqCountMatrix,annResult[,1:2])
htseqFPKMMatrixx.immune = createMatrixForImmuneInfiltration(htseqFPKMMatrix,annResult[,1:2])
htseqTPMMatrix.immune = createMatrixForImmuneInfiltration(htseqTPMMatrix,annResult[,1:2])


saveTSV(htseqCountMatrix.immune,output.counts.raw.infiltrationMatrixFileName)
saveTSV(htseqFPKMMatrixx.immune,output.counts.fpkm.infiltrationMatrixFileName)
saveTSV(htseqTPMMatrix.immune,output.counts.tpm.infiltrationMatrixFileName)



