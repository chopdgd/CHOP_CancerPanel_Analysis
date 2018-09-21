################################
# Author: Komal Rathi
# Purpose: Fix V8 minor changes
# Date: 02/16/2018
################################

library(reshape2)

# Read in cleaned data
data <- read.delim("../data/CleanDataFinal_V8.txt", stringsAsFactors = F)
colnames(data)[29] <- "Age_Years"
data$Patient_ID <- gsub('[a-zA-Z]$','',data$Study_ID)
data[which(data$Study_ID == "B5"),'Relapse_recurrent_tumor_specimen'] <- 'No'
data[which(data$Study_ID == "B9"),'Relapse_recurrent_tumor_specimen'] <- 'Yes'
data$CNV_Type[data$CNV_Type == "loss"] <- 'Loss'
data$CNV_Type[data$CNV_Type == "gain"] <- 'Gain'

# remove TPMT and NUDT15 from brain and solid tumors
data <- data[-which(data$Gene %in% c("TPMT","NUDT15") & data$Type %in% c("Brain Tumor","Solid Tumor (Non-CNS)")),]

# fix some variant tiers
data$Variant_Tier[data$Variant_Tier == "Tier 1B"] <- "1B"
data$Variant_Tier[data$Variant_Tier == "CNV Tier 1"] <- "CNV Tier 1A"
data$Variant_Tier[data$Variant_Tier == "Fusion Tier 1"] <- "Fusion Tier 1A"

write.table(data, "../data/CleanDataFinal_V9.txt", sep = "\t", row.names = F)
