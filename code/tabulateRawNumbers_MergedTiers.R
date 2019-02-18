################################
# Author: Komal Rathi
# Purpose: Tabulate raw numbers
# Date: 10/25/2018
# repeat the same but merging tier 1 and tier 2
# for the bigTable, combine Tier 1 and Tier 2 
################################

library(reshape2)
library(dplyr)

# Read in cleaned data
data <- read.delim("../data/CleanDataFinal_V10.txt")

# add tier type
data$Tier_Type <- ifelse(data$Variant_Tier %in% c("3","4","CNV Tier 3","Fusion Tier 3"), "Low", "High")
data$Tier_Type[data$Variant_Tier == "No Fusions Detected"] <- "NA"

# remove unnecessary
# remove TPMT/NUDT15/MPL from High Tiers CNS/Solid Tumors
# remove TPMT/NUDT15 from High Tiers Liquid Tumors
data <- data[-which(data$Type %in% c("Solid Tumor (Non-CNS)", "Brain Tumor") & data$Gene %in% c("TPMT","NUDT15","MPL") & data$Tier_Type == "High"),]
data <- data[-which(data$Type == "Liquid Tumor" & data$Gene %in% c("TPMT","NUDT15") & data$Tier_Type == "High"),]

np <- length(unique(data$Patient_ID))
ns <- length(unique(data$Sample_ID))

# combine Tiers
data$Variant_Tier <- as.character(data$Variant_Tier)
data$Variant_Tier[data$Variant_Tier %in% c("1A","1B","2")] <- "Combined_SNV_Tier_1A-1B-2"
data$Variant_Tier[data$Variant_Tier %in% c("Fusion Tier 1A", "Fusion Tier 2")] <- "Combined_Fusion_Tier_1A-2"
data$Variant_Tier[data$Variant_Tier %in% c("CNV Tier 1A","CNV Tier 1B","CNV Tier 2")] <- "Combined_CNV_Tier_1A-1B-2"

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
data %>% filter(Tier_Type == "High") %>% 
  dplyr::select(Sample_ID, Type, Age_Years) %>% 
  unique %>% group_by(Type) %>% 
  summarise(nsamples = n())

# low tiers
data %>% filter(Tier_Type == "Low") %>% 
  dplyr::select(Sample_ID, Type, Age_Years) %>% 
  unique %>% group_by(Type) %>% 
  summarise(nsamples = n())
#### samples level ####

#### patient level ####
# get age summaries from here for sheet1
tmpDat <- unique(data[,c("Patient_ID", "Type","Age_Years")])
tmpDat %>% group_by(Type) %>% 
  summarise(nsamples = n_distinct(Patient_ID), minAge = min(Age_Years), 
            maxAge = max(Age_Years), meanAge = mean(Age_Years),
            medianAge = median(Age_Years)) %>%
  mutate(age = paste0(minAge,'-',maxAge)) %>% 
  dplyr::select(-minAge, -maxAge)

# high tiers 
dataHT <- data[data$Tier_Type == "High",]
dataHT <- unique(dataHT[,c("Patient_ID", "Type","Age_Years")])
dataHT %>% group_by(Type) %>% summarise(nsamples = n_distinct(Patient_ID))

# low tiers
dataLT <- data[data$Tier_Type == "Low",]
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
write.table(res, file = '../tables/summary_sheet_per_tumortype_combined.txt', quote = F, sep = "\t", row.names = F)
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
write.table(res, file = '../tables/summary_sheet_all_tumortype_combined.txt', quote = F, sep = "\t", row.names = F)

# number of unique genes (n = 101 for Tab7, C24)
test <- data %>% group_by(Variant_Tier, Tier_Type) %>%
  filter(Gene != "") %>%
  summarise(ct = n(), genes = n_distinct(Gene))
test <- data %>%
  filter(Gene != "" & Variant_Tier == "Combined_SNV_Tier_1A-1B-2") %>%
  dplyr::select(Gene)
length(unique(test$Gene))

#### combining all ####

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

# now combine Tier1 and Tier2
tmp$Diagnostic_impact <- apply(tmp, 1, FUN = function(x) any(x[c('Diagnostic_impact_Tier_1', 'Diagnostic_impact_Tier_2')] == "Yes"))
tmp$Prognostic_impact <- apply(tmp, 1, FUN = function(x) any(x[c('Prognostic_impact_Tier_1', 'Prognostic_impact_Tier_2')] == "Yes"))
tmp$Diagnostic_impact <- ifelse(tmp$Diagnostic_impact == "FALSE", "No", "Yes")
tmp$Prognostic_impact <- ifelse(tmp$Prognostic_impact == "FALSE", "No", "Yes")

tmp$Diag_Prog_combined <- apply(tmp, 1, FUN = function(x) any(x[c('Diagnostic_impact', 'Prognostic_impact')] == "Yes"))
tmp$Diag_Prog_combined <- ifelse(tmp$Diag_Prog_combined == FALSE, "No", "Yes")

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
res <- res[,c("Type","Diagnostic_impact","Diagnostic_impact_Tier_1","Diagnostic_impact_Tier_2",
              "Therapeutic_impact_changed_diagnosis","Prognostic_impact","Prognostic_impact_Tier_1","Prognostic_impact_Tier_2",
              "Diag_Prog_combined",
              "Therapeutic_impact_based_on_risk","Therapeutic_impact_targeted","Suspected_germline_V")]

res.total <- countit(tmp)
res.total <- rbind(res.total, c("367","Total"))
rownames(res.total)[nrow(res.total)] <- "Type"
res.total$V2 <- NULL
res.total <- t(res.total)
res <- rbind(res, res.total)
write.table(res, file = '../tables/FinalBigTable_combined.txt', quote = F, sep = "\t", row.names = F)
#### big table ####

#### row totals for big table ####
tmp.row <- tmp
tmp.row$Suspected_germline_V <- NULL
tmp.row$results <- apply(tmp.row, 1, FUN = function(x) any(x[3:9] == "Yes"))
tmp.row <- tmp.row[,c(1,2,ncol(tmp.row))]
plyr::count(tmp.row,c('Type','results'))

################################
# TPMT and NUDT15 patient
tn <- data[which(data$Gene %in% c("TPMT","NUDT15")),]
p <- as.character(unique(tn$Patient_ID))
write.table(p, file = '../data/TPMT_NUDT15_Patients.txt', quote = F, row.names = F, col.names = F)
