########################################
# Author: Komal Rathi
# Purpose: Clean up V3 to get V4
# Date: 01/09/18
########################################

data <- read.delim("../data/CleanDataFinal_V3.txt"); # 4969 -> 5118 rows, 425 unique ids

# format Type column
data[,"Type"] <- gsub("NonCNS_Solid_v4_LFS", "noncns_solid", data[,"Type"])
data[,"Type"] <- gsub("LiquidTumors_v5_SPM", "liquidtumor", data[,"Type"])
data[,"Type"] <- gsub("Heme_nonmalignant_SPM", "hemenm_spm", data[,"Type"])
data[,"Type"] <- gsub("BrainTumorV4_LFS", "braintumor", data[,"Type"])
data[,"Type"] <- gsub("noncns_solid", "Solid Tumor (Non-CNS)", data[,"Type"])
data[,"Type"] <- gsub("liquidtumor", "Liquid Tumor", data[,"Type"])
data[,"Type"] <- gsub("hemenm_spm", "Non-Malignant Heme", data[,"Type"])
data[,"Type"] <- gsub("braintumor", "Brain Tumor", data[,"Type"])

# change CNV_Type Amp to Amplified
data$CNV_Type <- as.character(data$CNV_Type)
data$CNV_Type[data$CNV_Type == "Amp"] <- "Amplified"

# replace NAs with empty values
data[is.na(data)] <- ''

# fix Target_Position
data <- data %>% mutate(Target_Position = ifelse(Type_varType == "snp", gsub('-.*','',Target_Position), Target_Position))
data$Target_Position <- gsub('^Chr','chr',data$Target_Position)

# convert to character
data$Comments <- as.character(data$Comments)
data$Chromosome <- as.character(data$Chromosome)

# fix subcategory
data$Subcategory <- as.character(data$Subcategory)
data$Subcategory[data$Subcategory == "Lymphoma?"] <- 'Lymphoma'
data$Subcategory[data$Subcategory == "Embryonal \\blastomas\\\\\\"] <- 'Embryonal Blastomas'
data$Subcategory[data$Subcategory == "Germ Cell tumors"] <- 'Germ cell tumors'

# only for non empty Type_varType, mutate columns
data.uniq <- data %>% group_by(Sample_ID, Type, X3_Transcript, X5_Transcript, 
                               Amino_Acid_Change, CNV_Type, Chromosome_Band, 
                               Effect_codingEffect, Gene, HGVS_cDNA_level_nomenclature, 
                               HGVS_protein_level_nomenclature_pNomen3, Protein_Change, 
                               Start, Stop, Target_Position, Transcript, Type_varType, 
                               Variant_Tier, Study_ID, 
                               Type_of_Panel, Age__years_, Final_clinical_diagnosis, 
                               Overarching_Category, Subcategory, 
                               Primary_tumor_specimen, Relapse_recurrent_tumor_specimen, 
                               Therapeutic_impact_changed_diagnosis,
                               Summary_of_Findings, Diagnostic_impact_Tier_1, 
                               Diagnostic_impact_Tier_2,
                               Prognostic_impact_Tier_1, Prognostic_impact_Tier_2, 
                               Therapeutic_impact_targeted,
                               Therapeutic_impact_based_on_risk, Suspected_germline_V, 
                               Confirmed_germline_V) %>% 
  mutate(X_Mutation_in_Tumor = max(X_Mutation_in_Tumor), 
         CNV_Size = ifelse(Type_varType != "", toString(unique(CNV_Size)), CNV_Size),
         Chromosome = ifelse(Type_varType != "", toString(unique(Chromosome)), Chromosome),
         Comments = ifelse(Type_varType != "", toString(unique(Comments)), Comments)) %>%
  unique() %>% as.data.frame()

data.uniq$CNV_Size <- gsub('^, |, $', '', data.uniq$CNV_Size)
data.uniq$Comments <- gsub('^, |, $', '', data.uniq$Comments)
data.uniq$Chromosome <- gsub('^, |, $', '', data.uniq$Chromosome)

# remove chr from chromosome for consistency
data.uniq$Chromosome <- gsub('^chr','',data.uniq$Chromosome)

# empty values - sent to Lea
tmp <- data.uniq[which(data.uniq$Variant_Tier == ""),]
# write.table(tmp, file = '../data/moreMissingData.txt', sep = "\t", row.names = F)

# remove all "do not include" from missing data - fixed by Lea
missing <- read.delim('../data/moreMissingData_LFS Edits.txt', check.names = F)
setdiff(colnames(data.uniq),colnames(missing))
missing <- missing[which(missing$Variant_Tier != "do not include"),]
missing[is.na(missing)] <- ''
data.uniq <- data.uniq[which(data.uniq$Variant_Tier != ""),]
data.uniq <- rbind(data.uniq, missing) # 4869 rows, 424 sample IDs, 7 undetected (2 Fusions, 5 SNVs/CNVs)

# filter by age > 26
data.age.filt <- data.uniq %>% filter(Age__years_ <= 26) # 4780, 414 sample IDs, 4 undetected (2 each)

# write out version 3
write.table(data.age.filt, "../data/CleanDataFinal_V4.txt", sep = "\t", row.names = F)

