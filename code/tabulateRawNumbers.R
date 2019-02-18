################################
# Author: Komal Rathi
# Purpose: Tabulate raw numbers
# Date: 11/01/2018
################################

library(reshape2)
library(dplyr)
library(tidyr)

# Read in cleaned data
data <- read.delim("../data/CleanDataFinal_V10.txt", stringsAsFactors = F)

# add tier type
data$Tier_Type <- ifelse(data$Variant_Tier %in% c("3","4","CNV Tier 3","Fusion Tier 3"), "Low", "High")
data$Tier_Type[data$Variant_Tier == "No Fusions Detected"] <- "NA"

# remove unnecessary
# remove TPMT/NUDT15/MPL from High Tiers CNS/Solid Tumors
# remove TPMT/NUDT15 from High Tiers Liquid Tumors
data <- data[-which(data$Type %in% c("Solid Tumor (Non-CNS)", "Brain Tumor") & data$Gene %in% c("TPMT","NUDT15","MPL") & data$Tier_Type == "High"),]
data <- data[-which(data$Type == "Liquid Tumor" & data$Gene %in% c("TPMT","NUDT15") & data$Tier_Type == "High"),]

np <- length(unique(data$Patient_ID)) # 367
ns <- length(unique(data$Sample_ID)) # 389

# ave age across dataset
mean(data$Age_Years) # 8.580821
median(data$Age_Years) # 7

# patient IDs that have multiple sample IDs
ct <- unique(data[,c('Sample_ID','Patient_ID')])
ct <- plyr::count(ct$Patient_ID)
ct <- ct[which(ct$freq > 1),]
write.table(ct$x, file = '../data/PatientID_withMultipleSampIDs.txt', quote = F, sep = "\t", row.names = F, col.names = F)


#### samples level ####
# sheet1
tmpDat <- unique(data[,c("Sample_ID", "Type","Age_Years")])
tmpDat %>% group_by(Type) %>% 
  summarise(nsamples = n(), minAge = min(Age_Years), 
            maxAge = max(Age_Years), meanAge = mean(Age_Years),
            medianAge = median(Age_Years)) %>%
  mutate(age = paste0(minAge,'-',maxAge)) %>% 
  dplyr::select(-minAge, -maxAge)

# high tiers 
highTiers <- c("1A", "1B", "2", 
               "CNV Tier 1A", "CNV Tier 1B", "CNV Tier 2", 
               "Fusion Tier 1A", "Fusion Tier 2")
dataHT <- data[data$Variant_Tier %in% highTiers,]
dataHT <- unique(dataHT[,c("Sample_ID", "Type","Age_Years")])
dataHT %>% group_by(Type) %>% summarise(nsamples = n())

# low tiers
lowTiers <- c("3", "4", "CNV Tier 3", "Fusion Tier 3") 
dataLT <- data[data[,"Variant_Tier"]%in%lowTiers,]
dataLT <- unique(dataLT[,c("Sample_ID", "Type","Age_Years")])
dataLT %>% group_by(Type) %>% summarise(nsamples = n())
#### samples level ####

#### patient level ####
# get age summaries from here for sheet1
tmpDat <- unique(data[,c("Patient_ID", "Type","Age_Years")])
tmpDat %>% group_by(Type) %>% 
  summarise(nsamples = n_distinct(Patient_ID), minAge = min(Age_Years), 
            maxAge = max(Age_Years), meanAge = mean(Age_Years),
            medianAge = median(Age_Years)) %>%
  mutate(age = paste0(minAge,'-',maxAge)) %>% dplyr::select(-minAge, -maxAge)

# high tiers 
highTiers <- c("1A", "1B", "2", 
               "CNV Tier 1A", "CNV Tier 1B", "CNV Tier 2", 
               "Fusion Tier 1A", "Fusion Tier 2")
dataHT <- data[data$Variant_Tier %in% highTiers,]
dataHT <- unique(dataHT[,c("Patient_ID", "Type","Age_Years")])
dataHT %>% group_by(Type) %>% summarise(nsamples = n_distinct(Patient_ID))

# low tiers
lowTiers <- c("3", "4", "CNV Tier 3", "Fusion Tier 3") 
dataLT <- data[data[,"Variant_Tier"]%in%lowTiers,]
dataLT <- unique(dataLT[,c("Patient_ID", "Type","Age_Years")])
dataLT %>% group_by(Type) %>% summarise(nsamples = n_distinct(Patient_ID))
#### patient level ####

#### per Type per Tier ####
tmp <- unique(data[,c('Patient_ID','Sample_ID','Type','Tier_Type','Variant_Tier','Age_Years')])

# ages of the samples
with.ages <- tmp %>% group_by(Type, Variant_Tier) %>% 
  summarise(nsamples = n(), 
            minAge = min(Age_Years), maxAge = max(Age_Years), 
            Mean_Age = round(mean(Age_Years),2), Median_Age = median(Age_Years), 
            nPatients = n_distinct(Patient_ID)) %>% 
  mutate(age = paste0(minAge,'-',maxAge), 
         totalPatients = sum(nPatients),
         perc = round(nPatients/totalPatients*100,2)) %>% 
  dplyr::select(-minAge, -maxAge)

# add variant numbers
vars <- data %>% group_by(Type, Variant_Tier, Tier_Type) %>% summarise(ct = n())
res <- merge(with.ages, vars, by = c("Type",'Variant_Tier'))
res <- res[,c('Type','Variant_Tier','Tier_Type','ct','nsamples','nPatients','perc','age','Mean_Age','Median_Age')]
colnames(res) <- c("Type","Variant_Tier","Tier_Type","# of Variants","# of unique samples per Variant Tier","Unique Patient ID Per Tier","Patient per Tier (%)","Age_range","Mean_Age","Median_Age")
write.table(res, file = '../tables/summary_sheet_per_tumortype.txt', quote = F, sep = "\t", row.names = F)
#### per Type per Tier ####

#### combining all ####
tmp <- unique(data[,c('Patient_ID','Sample_ID','Variant_Tier','Age_Years')])
with.ages <- tmp %>% group_by(Variant_Tier) %>% 
  summarise(nsamples = n_distinct(Sample_ID), 
            minAge = min(Age_Years), maxAge = max(Age_Years), 
            Mean_Age = round(mean(Age_Years),2), Median_Age = median(Age_Years), 
            nPatients = n_distinct(Patient_ID)) %>% 
  mutate(age = paste0(minAge,'-',maxAge), 
         totalPatients = sum(nPatients),
         perc = round(nPatients/totalPatients*100,2)) %>% 
  dplyr::select(-minAge, -maxAge)

vars <- data %>% group_by(Variant_Tier, Tier_Type) %>% summarise(ct = n())
res <- merge(with.ages, vars, by = c('Variant_Tier'))
res <- res[,c('Variant_Tier','Tier_Type','ct','nsamples','nPatients','perc','age','Mean_Age','Median_Age')]
colnames(res) <- c("Variant_Tier","Tier_Type","# of Variants","# of unique samples per Variant Tier","Unique Patient ID Per Tier","Patient per Tier (%)","Age_range","Mean_Age","Median_Age")
write.table(res, file = '../tables/summary_sheet_all_tumortype.txt', quote = F, sep = "\t", row.names = F)
#### combining all ####

#### Number of primary and relapse samples ####
tmp <- unique(data[,c('Sample_ID','Patient_ID','Primary_tumor_specimen','Relapse_recurrent_tumor_specimen')])
tmp$Primary_tumor_specimen <- as.character(tmp$Primary_tumor_specimen)
tmp$Relapse_recurrent_tumor_specimen <- as.character(tmp$Relapse_recurrent_tumor_specimen)
tmp$Primary_tumor_specimen <- gsub('Yes.*|ACUTE.*','Yes',tmp$Primary_tumor_specimen)
tmp$Relapse_recurrent_tumor_specimen <- gsub('Yes.*','Yes',tmp$Relapse_recurrent_tumor_specimen)
tmp$Relapse_recurrent_tumor_specimen <- gsub("No.*","No", tmp$Relapse_recurrent_tumor_specimen)

tmpp <- unique(tmp[,c('Sample_ID','Patient_ID','Primary_tumor_specimen','Relapse_recurrent_tumor_specimen')])
table(tmpp$Primary_tumor_specimen)
table(tmpp$Relapse_recurrent_tumor_specimen)
#### Number of primary and relapse samples ####

#### big table ####
tmp <- data[,c("Patient_ID","Type",
               "Diagnostic_impact_Tier_1","Diagnostic_impact_Tier_2",
               "Therapeutic_impact_changed_diagnosis",
               "Prognostic_impact_Tier_1","Prognostic_impact_Tier_2",
               "Therapeutic_impact_based_on_risk","Therapeutic_impact_targeted",
               "Suspected_germline_V")]
tmp <- unique(tmp)

# for number of  total cases
plyr::count(unique(tmp[,c("Patient_ID","Type")]),"Type")

fixit <- function(x){
  x <- gsub("Yes.*|WHSCH1.*|PAX5.*|CDKN2A.*","Yes", x)
  x <- gsub("No.*","No", x)
  return(x)
}

tmp[,3:ncol(tmp)] <- apply(tmp[,3:ncol(tmp)], MARGIN = 2, FUN = function(x) fixit(x))
apply(tmp[,3:ncol(tmp)], MARGIN = 2, FUN = function(x) unique(x))

countit <- function(x){
  len <- length(unique(x$Patient_ID))
  y <- apply(x[,3:ncol(x)], 2, FUN = function(x) length(unique(grep('Yes',x))))
  y <- as.data.frame(y)
  colnames(y) <- 'V1'
  y$V2 <- round(y$V1/len*100,2)
  y$V1 <- paste0(y$V1,' (',y$V2,'%)')
  y$V2 <- rownames(y)
  return(y)
}

res <- plyr::ddply(tmp, .variables = 'Type', .fun = function(x) countit(x))
res <- dcast(res, Type~V2, value.var = "V1")
res <- res[,c("Type","Diagnostic_impact_Tier_1","Diagnostic_impact_Tier_2",
              "Therapeutic_impact_changed_diagnosis","Prognostic_impact_Tier_1","Prognostic_impact_Tier_2",
              "Therapeutic_impact_based_on_risk","Therapeutic_impact_targeted","Suspected_germline_V")]

res.total <- countit(tmp)
res.total <- rbind(res.total, c("367","Total"))
rownames(res.total)[nrow(res.total)] <- "Type"
res.total$V2 <- NULL
res.total <- t(res.total)
res <- rbind(res, res.total)
write.table(res, file = '../tables/FinalBigTable_tab1.txt', quote = F, sep = "\t", row.names = F)
#### big table ####

#### row totals for big table ####
tmp.row <- tmp
tmp.row$Suspected_germline_V <- NULL
tmp.row$results <- apply(tmp.row, 1, FUN = function(x) any(x[3:9] == "Yes"))
tmp.row <- tmp.row[,c(1,2,ncol(tmp.row))]
plyr::count(tmp.row,c('Type','results'))

