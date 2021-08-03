#TODO de donde sale MSI status y los tratamientos


htseq = HELPER_LoadTSV("./resources/human/gdc/tcga/coad/rawData/HTSeq-FPKM-Samples.tsv")
# remove with 2 replicas
htseq.duplicated = htseq$`Sample ID`[duplicated(htseq$`Sample ID`)]
htseq.final = subset(htseq, !htseq$`Sample ID` %in% htseq.duplicated)
#only keep primary tumor and control
htseq.final = htseq.final[htseq.final$`Sample Type` == "Primary Tumor" | htseq.final$`Sample Type` == "Solid Tissue Normal",]
colnames(htseq.final)[8] = "GROUP"
htseq.final = htseq.final[order(htseq.final$GROUP, htseq.final$`Sample ID`),]

#sacamos de primary tumor los que tienen mas de 1
onlyTumor = htseq.final[htseq.final$GROUP == "Primary Tumor",]
repetead = onlyTumor[duplicated(onlyTumor$`Case ID`),]

htseq.final= subset(htseq.final, !htseq.final$`Sample ID` %in% repetead$`Sample ID`)

table(htseq.final$GROUP) #496 final = 456 primary tumor y 41 solid tissue normal
table(duplicated(htseq.final$`Sample ID`)) #0 duplicados
colnames(htseq.final) = gsub(" ","_", colnames(htseq.final))
HELPER_SaveCSV(htseq.final, "./resources/human/gdc/tcga/coad/TCGA_COAD_HTSeq-FPKM_samples.csv", rowNames = FALSE)

#********


clinical = HELPER_LoadTSV("./resources/human/gdc/tcga/coad/rawData/clinical/clinical.tsv")
exposure = HELPER_LoadTSV("./resources/human/gdc/tcga/coad/rawData/clinical/exposure.tsv")
all_CDEs = HELPER_LoadTSV("./resources/human/gdc/tcga/coad/rawData/clinical/all_cdes.txt")
auxiliary = HELPER_LoadTSV("./resources/human/gdc/tcga/coad/rawData/clinical/auxiliary_public_coad.txt")
clinical_drug = HELPER_LoadTSV("./resources/human/gdc/tcga/coad/rawData/clinical/clinical_drug_public_coad.txt")
coadr.merged = HELPER_LoadTSV("./resources/human/gdc/tcga/coad/rawData/clinical/COAD.merged_only_auxiliary_clin_format.txt")

exposure = exposure[,-(2:3)]

clinical.final = merge(clinical, exposure, by.x = "case_id", by.y = "case_id" )


msi_status = as.data.frame(t(coadr.merged))
msi_status = msi_status[-1,]
msi_status = data.frame(case_id = as.character(msi_status$V11), msi_status = as.character(msi_status$V13), stringsAsFactors = FALSE)

clinical.final = merge(clinical.final, msi_status, by.x = "case_id", by.y = "case_id", all.x = TRUE )



plus = as.data.frame(t(all_CDEs))
plus = plus[-1,]


plus = data.frame (case_id = plus$V1, anatomic_neoplasm_subdivision = plus$V8, braf_gene_analysis_result = plus$V19, colon_polyps_present =  plus$V26,histological_type =  plus$V44, history_of_colon_polyps = plus$V45, followup_treatment_success = plus$V40,
kras_gene_analysis_performed = plus$V53, kras_mutation_codon = plus$V54, kras_mutation_found = plus$V55, loss_expression_of_mismatch_repair_proteins_by_ihc = plus$V56, loss_expression_of_mismatch_repair_proteins_by_ihc_result = plus$V57,
loss_expression_of_mismatch_repair_proteins_by_ihc_result = plus$V58, loss_expression_of_mismatch_repair_proteins_by_ihc_result = plus$V59, loss_expression_of_mismatch_repair_proteins_by_ihc_result = plus$V60,
lymph_node_examined_count = plus$V62, lymphatic_invasion = plus$V63,
number_of_lymphnodes_positive_by_he = plus$V74, number_of_lymphnodes_positive_by_ihc = plus$V75 , number_of_first_degree_relatives_with_cancer_diagnosis = plus$V72, person_neoplasm_cancer_status = plus$V82,
postoperative_rx_tx = plus$V83,primary_therapy_outcome_success = plus$V87,
synchronous_colon_cancer_present = plus$V96, venous_invasion = plus$V103)
plus = plus[-1,]

clinical.final = merge(clinical.final, plus, by.x = "case_id", by.y = "case_id", all.x = TRUE )



clinical.final = clinical.final[clinical.final$submitter_id %in% htseq.final$Case_ID,]
dim(clinical.final)
clinical.final$year_of_birth = as.numeric(as.character(clinical.final$year_of_birth))
clinical.final$year_of_death = as.numeric(as.character(clinical.final$year_of_death))
clinical.final$age_at_diagnosis = as.numeric(as.character(clinical.final$age_at_diagnosis))
clinical.final$days_to_birth = as.numeric(as.character(clinical.final$days_to_birth))
clinical.final$days_to_death = as.numeric(as.character(clinical.final$days_to_death))
clinical.final$days_to_last_follow_up = as.numeric(as.character(clinical.final$days_to_last_follow_up))
clinical.final$days_to_last_known_disease_status = as.numeric(as.character(clinical.final$days_to_last_known_disease_status))
clinical.final$days_to_recurrence = as.numeric(as.character(clinical.final$days_to_recurrence))
clinical.final$cigarettes_per_day = as.numeric(as.character(clinical.final$cigarettes_per_day))
clinical.final$years_smoked = as.numeric(as.character(clinical.final$years_smoked))
clinical.final$bmi = as.numeric(as.character(clinical.final$bmi))
clinical.final$weight = as.numeric(as.character(clinical.final$weight))
clinical.final$height = as.numeric(as.character(clinical.final$height))



colnames(clinical.final) = gsub(" ","_", colnames(clinical.final))
colnames(clinical_drug) = gsub(" ","_", colnames(clinical_drug))

HELPER_SaveCSV(clinical.final, "./resources/human/gdc/tcga/coad/TCGA_COAD_ClinicalMerge.csv", rowNames = FALSE)
HELPER_SaveCSV(clinical_drug, "./resources/human/gdc/tcga/coad/TCGA_COAD_ClinicalDrugs.csv", rowNames = FALSE)

#*******
onlyTumor = htseq.final[htseq.final$GROUP == "Primary Tumor",]
colnames(onlyTumor)[6] = "Case_ID"
#agrego las columnas de clinical para poder separar por grupo
onlyTumor.final = merge(onlyTumor, clinical.final, by.x = "Case_ID" , by.y = "submitter_id")
colnames(onlyTumor.final) = gsub(" ","_", colnames(onlyTumor.final))

HELPER_SaveCSV(onlyTumor.final, "./resources/human/gdc/tcga/coad/TCGA_COAD_Only_tumor_HTSeq-FPKM_samples.csv", rowNames = FALSE)
View(onlyTumor.final)




