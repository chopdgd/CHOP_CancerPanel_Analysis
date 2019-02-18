########################################
#Author: Komal Rathi
#Purpose: Create Supplemental Table
#Date: 10/29/2018
########################################

# read V10 and subset to a specific set of columns
dat <- read.delim('../data/CleanDataFinal_V10.txt', stringsAsFactors = F)
cols <- c('Study_ID','Age_Years','Type','Subcategory','Primary_tumor_specimen',
          'Gene','Transcript','Amino_Acid_Change','Protein_Change','X_Mutation_in_Tumor',
          'Chromosome','Chromosome_Band','CNV_Type','Variant_Tier',
          'Diagnostic_impact_Tier_1','Diagnostic_impact_Tier_2','Therapeutic_impact_changed_diagnosis',
          'Prognostic_impact_Tier_1','Prognostic_impact_Tier_2','Therapeutic_impact_based_on_risk',
          'Therapeutic_impact_targeted','Suspected_germline_V')
dat <- dat[,cols]

# replace samples where only Tier 4 was detected with No variants detected
x <- plyr::count(dat, c("Study_ID"))
x <- x[which(x$freq == 1),]
s.to.replace <- dat[which(dat$Study_ID %in% x$Study_ID),]
s.to.replace <- s.to.replace[which(s.to.replace$Variant_Tier == "4"),'Study_ID']
dat[which(dat$Study_ID %in% s.to.replace),'Variant_Tier'] <- 'No Variants Detected'

# remove all Tier 4 variants where multiple variant types were detected
x <- unique(dat[which(dat$Variant_Tier %in% "4"),'Study_ID'])
x <- dat[which(dat$Study_ID %in% x),]
x <- plyr::count(x$Study_ID)
if(nrow(x[which(x$freq == 1),]) == 0){
  print("None found")
  dat <- dat[which(dat$Variant_Tier != "4"),]
} # this should be 0

# Now merge Prognostic/Diagnostic Tier 1 and Tier 2 columns
# To do that, first remove all extra information in the columns
fixit <- function(x){
  x <- gsub("Yes.*|WHSCH1.*|PAX5.*|CDKN2A.*","Yes", x)
  x <- gsub("No.*","No", x)
  return(x)
}
dat[,15:ncol(dat)] <- apply(dat[,15:ncol(dat)], MARGIN = 2, FUN = function(x) fixit(x))
apply(dat[,15:ncol(dat)], MARGIN = 2, FUN = function(x) unique(x))
dat$Diagnostic_impact <- apply(dat, 1, FUN = function(x) any(x[c('Diagnostic_impact_Tier_1', 'Diagnostic_impact_Tier_2')] == "Yes"))
dat$Prognostic_impact <- apply(dat, 1, FUN = function(x) any(x[c('Prognostic_impact_Tier_1', 'Prognostic_impact_Tier_2')] == "Yes"))
dat$Diagnostic_impact <- ifelse(dat$Diagnostic_impact == "FALSE", "No", "Yes")
dat$Prognostic_impact <- ifelse(dat$Prognostic_impact == "FALSE", "No", "Yes")

cols <- c('Study_ID','Age_Years','Type','Subcategory','Primary_tumor_specimen',
          'Gene','Transcript','Amino_Acid_Change','Protein_Change','X_Mutation_in_Tumor',
          'Chromosome','Chromosome_Band','CNV_Type','Variant_Tier',
          'Diagnostic_impact','Therapeutic_impact_changed_diagnosis',
          'Prognostic_impact','Therapeutic_impact_based_on_risk',
          'Therapeutic_impact_targeted','Suspected_germline_V')
dat <- dat[,cols]
dat[dat == ""] <- NA
apply(dat, MARGIN = 2, FUN = function(x) unique(x))
new.cols <- c("Study ID","Age, Years","Tumor Category","Tumor Subtype","Primary Tumor","Gene","Transcript","Coding change","Protein change","Variant Allele Frequency","Chromosome","CNV arm","CNV Change","Tier","Tier 1 or 2 Diagnostic Impact","Changed Diagnosis","Prognostic Impact","Based on Prognosis","Targeted Therapy Availability","Potential Germline Alteration")
colnames(dat) <- new.cols
write.table(dat, file = '../data/CleanDataFinal_V10_Supplemental.txt',  sep = "\t", row.names = F)
