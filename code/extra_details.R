########################################
#Author: Komal Rathi
#Purpose: Extra details for Lea/Suzanne
#Date: 10/29/2018
########################################

library(reshape2)
library(dplyr)

# Read in cleaned data
data <- read.delim("../data/CleanDataFinal_V10.txt", stringsAsFactors = F)

# add tier type
data$Tier_Type <- ifelse(data$Variant_Tier %in% c("3","4","CNV Tier 3","Fusion Tier 3"), "Low", "High")
data$Tier_Type[data$Variant_Tier == "No Fusions Detected"] <- "NA"

# how many liquid patients (not panels/tests) 
# had Tier 1B TPMT or NUDT15 variants
samples <- data %>% filter(Type == "Liquid Tumor" & 
                             Tier_Type == "High" & 
                             Gene %in% c('NUDT15','TPMT')) %>% 
  select(Patient_ID) %>% 
  unique()
write.table(samples, file = '../tables/LiquidTumors_NUDT15-TPMT_HighTiers.txt', quote = F, sep = "\t", row.names = F)

# remove some
data <- data[-which(data$Type %in% c("Solid Tumor (Non-CNS)", "Brain Tumor") & data$Gene %in% c("TPMT","NUDT15","MPL") & data$Tier_Type == "High"),]
data <- data[-which(data$Type == "Liquid Tumor" & data$Gene %in% c("TPMT","NUDT15") & data$Tier_Type == "High"),]
data <- unique(data)

# count the number of low tier specific panels out of 389
tmp <- data[,c('Patient_ID','Sample_ID','Tier_Type'),]
tmp <- plyr::count(tmp[,c('Sample_ID','Tier_Type')])
tmp[which(tmp$Tier_Type == 'NA'),2] <- 'Not_Avaiable'
tmp <- dcast(tmp, Sample_ID ~ Tier_Type, value.var = 'freq')
tmp[is.na(tmp)] <- 0
tmp <- tmp[which(tmp$High == 0 & tmp$Not_Avaiable == 0),]
samples <- unique(tmp$Sample_ID)
unique(data[which(data$Sample_ID %in% samples),'Tier_Type'])
write.table(samples, file = '../tables/Samples_with_LowTiers_only.txt', quote = F, sep = "\t", row.names = F)

# stats for the table in Raw Data Table
Input =("
Tumor.Type Present Not.Present
CNS 87 9
Leukemia 110 14
Solid 92 55")
Matriz = as.matrix(read.table(textConnection(Input),
                              header=TRUE, 
                              row.names=1))
chisq.test(Matriz) # p-value 4.612e-09
chisq.test(Matriz[c(1:2),]) # CNS vs Leukemia = 0.8116
chisq.test(Matriz[c(2,3),]) # Leukemia vs Solid = 1.768e-06
chisq.test(Matriz[c(1,3),]) # CNS vs Solid = 2.572e-06


Input =("
Tumor.Type Present Not.Present
CNS 5 91
Leukemia 7 117
Solid 32 115")
Matriz = as.matrix(read.table(textConnection(Input),
                              header=TRUE, 
                              row.names=1))
chisq.test(Matriz) # p-value 1.484e-05
chisq.test(Matriz[1:2,]) # CNS vs Leukemia = 1
chisq.test(Matriz[2:3,]) # Leukemia vs Solid = 0.000326
chisq.test(Matriz[c(1,3),]) # Solid vs CNS = 0.0008684
