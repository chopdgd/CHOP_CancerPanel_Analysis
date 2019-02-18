########################################
#Author: Pichai Raman
#Purpose: Clean data
#Date: 10/02/17
########################################

#Call libraries
library("tidyverse")
data <- read.delim("../data/mergedDataFromRaw.txt");

#Clean up column names
colnames(data) <- gsub("\\.", "_", colnames(data))
colnames(data) <- gsub("__", "_", colnames(data))
colnames(data)[c(48, 120)]<- c("Effect_codingEffect", "Type_varType")

#Read Marilyn updates and merge with data
updateData <- read.delim("../data/ReportedFromFM.txt");
colnames(updateData) <- gsub("\\.", "_", colnames(updateData))
colnames(updateData) <- gsub("__", "_", colnames(updateData))
colnames(updateData)[c(12, 21)]<- c("Effect_codingEffect", "Type_varType")
updateData <- updateData[1:22];
mapping <- read.delim("../data/mapping2.txt")
colnames(mapping) <- c("Type", "Sample_ID")

mapping[,"Sample_ID"] <- gsub(" ", "", mapping[,"Sample_ID"])
updateData[,"Sample_ID"] <- gsub(" ", "", updateData[,"Sample_ID"])
updateData <- merge(mapping, updateData, by.x="Sample_ID", by.y="Sample_ID")
updateData <- updateData[,c(2, 1, 4:23)]
colnames(updateData)[1] <- "Type"

data <- rbind(data[,colnames(updateData)], updateData)

#Add in last piece of data / disease
brainData <- read.delim("../data/manuscript_samples_data_combined v4_0811_DEIDENTIFIED_Brain.txt")
hemeData <- read.delim("../data/manuscript_samples_data_combined v4_0811_DEIDENTIFIED_HemeNonMalignant.txt")
liquidData <- read.delim("../data/manuscript_samples_data_combined v4_0811_DEIDENTIFIED_LiquidTUmor.txt")
solidData <- read.delim("../data/manuscript_samples_data_combined v4_0811_DEIDENTIFIED_solid.txt")

pullCols <- function(x)
{
	keepCols <- c("Order.Number", "Study.ID", "Type.of.Panel", "Age..years.", "Final.clinical.diagnosis", "Overarching.Category", "Subcategory", "Primary.tumor.specimen", "Relapse.recurrent.tumor.specimen", "Therapeutic.impact.changed.diagnosis");
	x <- x[,keepCols]
	colnames(x) <- gsub("\\.", "_", colnames(x));
	colnames(x) <- gsub(" ", "_", colnames(x));
	colnames(x)[4] <- "Age_Years"
	x <- unique(x);
	return(x);
}
brainData <- pullCols(brainData);
hemeData <- pullCols(hemeData);
liquidData <- pullCols(liquidData);
solidData <- pullCols(solidData);
allDataDiag <- rbind(brainData, hemeData, liquidData, solidData)
data <- merge(data, allDataDiag, by.x="Sample_ID", by.y="Order_Number", all.x=T)

#Columns for Marilyn
data[,"ValidateRow"] <- F
data[,"ReasonToValidateRow"] <- NA


#Clean up each column
# [1] "Type": 4 levels, Brain, CNS, Solid, Liquid...broad category                                  
# [2] "Sample_ID": 418 levels - Patient / Sample ID                             
# [3] "X_Mutation_in_Tumor": Numeric field between 0 & .49, 2700 NA's                
# [4] "X3_Transcript": 73 levels, 5501 NA's                          
# [5] "X5_Transcript":  76 levels, 5499 NA's                          
# [6] "Amino_Acid_Change": 2275 levels, 2528 NA's                   
# [7] "CNV_Size": Numeric field between -3335- & 162400000, 2991 NA's                               
# [8] "CNV_Type": 51 levels, 3886 NA's : MUST CLEAN
data[,"CNV_Type"] <- as.character(data[,"CNV_Type"]);
data[,"CNV_Type"] <- toupper(data[,"CNV_Type"])
data[grepl("GAIN/AMP", data[,"CNV_Type"]), "CNV_Type"] <- "Gain/Amplified";
data[grepl("GAIN", data[,"CNV_Type"]), "CNV_Type"] <- "Gain";
data[grepl("LOSS", data[,"CNV_Type"]), "CNV_Type"] <- "Loss";
data[grepl("LOH", data[,"CNV_Type"]), "CNV_Type"] <- "LOH";
data[grepl("AMP", data[,"CNV_Type"]), "CNV_Type"] <- "Amplified"; 
data[grepl("COMPLEX CNV", data[,"CNV_Type"]), "CNV_Type"] <- "Complex CNV";
data[data[,"CNV_Type"]=="Gain/Amplified", "ValidateRow"] <- T;
data[data[,"CNV_Type"]=="Gain/Amplified", "ReasonToValidateRow"] <- "Gain/Amplified";
data[data[,"CNV_Type"]=="Gain/Amplified", "ValidateRow"] <- T;
data[data[,"CNV_Type"]=="Complex CNV", "ReasonToValidateRow"] <- "Complex CNV";

# [9] "Chormosome_Band" : MUST CLEAN -> 3 levels and 3561 NA's
colnames(data)[9] <- "Chromosome_Band"
data[,"Chromosome_Band"] <- toupper(data[,"Chromosome_Band"])
data[,"Chromosome_Band"] <- as.character(data[,"Chromosome_Band"]);
data[grepl("P", data[,"Chromosome_Band"]), "Chromosome_Band"] <- "p";
data[grepl("Q", data[,"Chromosome_Band"]), "Chromosome_Band"] <- "q";
data[grepl("WC", data[,"Chromosome_Band"]), "Chromosome_Band"] <- "Whole Chromosome";
data[grepl("WHOLE", data[,"Chromosome_Band"]), "Chromosome_Band"] <- "Whole Chromosome";
data[grepl("LONG", data[,"Chromosome_Band"]), "Chromosome_Band"] <- "q";

# [10] "Chromosome" : MUST CLEAN
data[,"Chromosome"] <- as.character(data[,"Chromosome"]);
data[,"Chromosome"] <- gsub("chr", "", data[,"Chromosome"])
data[,"Chromosome"] <- gsub(" ", "", data[,"Chromosome"])
data[,"Chromosome"] <- gsub("q", "", data[,"Chromosome"])
data[,"Chromosome"] <- gsub("p", "", data[,"Chromosome"])
data[,"Chromosome"] <- gsub("WC", "", data[,"Chromosome"])
data[,"Chromosome"] <- gsub(";", ",", data[,"Chromosome"])
data[,"Chromosome"] <- gsub("\\|", ",", data[,"Chromosome"])
data[grepl("7,10,16", data[,"Chromosome"]), "Chromosome"] <- "7,10,16";
data[grepl("1,2,3,7", data[,"Chromosome"]), "Chromosome"] <- "1,2,3,7,10,12,13,18";
data[data[,"Chromosome"]=="v", "ValidateRow"] <- T;
data[data[,"Chromosome"]=="v", "ReasonToValidateRow"] <- "Chromosome number v";

# [11] "Comments": Not Examined                              
# [12] "Effect_codingEffect": 2374 NA's, 2272 levels
# [13] "Gene": MUST CLEAN, a number of things don't have values
data[,"Gene"] <- as.character(data[,"Gene"]);

#Find anything with multiple genes and create more rows
data[grepl(",", data[,"Gene"]),"Gene"] <- "";
data[grepl(";", data[,"Gene"]), "Gene"] <- "";
data[grepl("\\|", data[,"Gene"]), "Gene"] <- "";
data[grepl(" of ", data[,"Gene"]), "Gene"] <- "";
data[grepl("no cancer ", data[,"Gene"]), "Gene"] <- "";
data[grepl("/", data[,"Gene"]), "Gene"] <- "";
data[grepl("HTATIP2 ST2", data[,"Gene"]),"Gene"]<- "";
data[grepl("HRAS CCND1", data[,"Gene"]), "Gene"] <- "";
data[grepl("FSTL3 RPS15 TCF3", data[,"Gene"]), "Gene"] <- "";
data[grepl("PTEN FAS", data[,"Gene"]), "Gene"] <- "";
data[grepl("PMS2 RAC1", data[,"Gene"]), "Gene"] <- "";
data[grepl("\\...", data[,"Gene"]), "Gene"] <- "";

data[grepl("\\)", data[,"Gene"]), "Gene"] <- "";
data[,"Gene"] <- gsub(" ", "", data[,"Gene"])
data[,"Gene"] <- gsub("\\.", "", data[,"Gene"])
data[,"Gene"] <- gsub(",", "", data[,"Gene"])
data[grepl("c-MYC", data[,"Gene"]), "Gene"] <- "MYC";

# [14] "HGVS_cDNA_level_nomenclature": 3979 NA's, 1342 levels : No change         
# [15] "HGVS_protein_level_nomenclature_pNomen3": 4055 NA's, 1232 levels : No change 
# [16] "Protein_Change": 2609 NA's, 2343 levels : No change                           
# [17] "Start": 4055 NA's, 1232 levels : No change                                   
# [18] "Stop": 4055 NA's, 1232 levels : No change                                   
# [19] "Target_Position" :882 NA's, 2814 levels : No change                         
# [20] "Transcript" :2491 NA's, 400 levels : No change                            
# [21] "Type_varType" : MUST CHANGE 
data[,"Type_varType"] <- as.character(data[,"Type_varType"])                         
data[grepl("ENSG", data[,"Type_varType"]), "Type_varType"] <- "";

#Convert delins --> indel
# [22] "Variant_Tier" 
data[,"Variant_Tier"] <- as.character(data[,"Variant_Tier"])   
data[grepl("not report", data[,"Variant_Tier"]), "Variant_Tier"] <- "do not report";
data[grepl("non-reportable", data[,"Variant_Tier"]), "Variant_Tier"] <- "do not report";
data[grepl("Fusion Tier 1A\n", data[,"Variant_Tier"]), "Variant_Tier"] <- "Fusion Tier 1A";

#Report anything that doesn't have a tier
mutTiers <- c("1A", "1B", "2", "3", "4")
data[data[,"Variant_Tier"]=="", "ValidateRow"] <- T;
data[data[,"Variant_Tier"]=="", "ReasonToValidateRow"] <- "No Variant Tier Listed";
data[data[,"Gene"]==""&data[,"Variant_Tier"]%in%mutTiers, "ValidateRow"] <- T;
data[data[,"Gene"]==""&data[,"Variant_Tier"]%in%mutTiers, "ReasonToValidateRow"] <- "Variant but no Gene Name";

#TMP dataMaster <- data; --> till I check that it all worked

#Let's try to filter out the rows that have "no variant tier listed"
#dim(data[data[,"ReasonToValidateRow"]%in%("No Variant Tier Listed"),]) # 417 rows currently
dataRemNoVariantTier <- data[!(data[,"ReasonToValidateRow"]%in%c("No Variant Tier Listed")), ]
dataNoVariantTier <- data[(data[,"ReasonToValidateRow"]%in%c("No Variant Tier Listed")), ]

#Remove ones that have overexpress
dataNoVariantTier <- dataNoVariantTier %>% filter(!(grepl("express", Comments))) # Now at 409 rows, removed 8 rows
dataNoVariantTier <- dataNoVariantTier %>% filter(!((!is.na(CNV_Size))&CNV_Size>1000)) # Now at 221 rows, removed 188 rows

dataMuts <- dataRemNoVariantTier %>% filter(Variant_Tier%in%c("1", "2", "3", "4"))
sampMutExist <- unique(paste(dataMuts[,"Sample_ID"], dataMuts[,"Amino_Acid_Change"], sep="_"));
dataNoVariantTier[,"SampMut"] <- paste(dataNoVariantTier[,"Sample_ID"], dataNoVariantTier[,"Amino_Acid_Change"], sep="_")
dataNoVariantTier <- dataNoVariantTier %>% filter(!(SampMut%in%sampMutExist)) # Now only 165 left, removed 56
dataNoVariantTier <- dataNoVariantTier[1:33]

data <- rbind(dataRemNoVariantTier, dataNoVariantTier);


write.table(data, "../data/CleanDataForPlotting.txt", sep="\t", row.names=F)


#Save this as RData
saveRDS(data, "../data/data.rds")
















