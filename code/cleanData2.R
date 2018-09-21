########################################
#Author: Pichai Raman
#Purpose: After a manual clean up with Marilyn and team, doing one last clean up
#Basically need to 1) parse out multichromosome values 2) remove any duplicates 
#Date: 10/02/17
########################################

#Call libraries
library("tidyverse")
data <- read.delim("../data/CleanDataForPlotting_update.txt");
dim(data) # 5091 / 35 

#First let's remove anything that says remove and get rid of cleaning columns
data <- data[data[,"Updated"]!="Remove",]
data <- data[1:31]; # 5060 / 31

#1 ) Now let's separate into multiple rows any that have multiple chromosomes
dataMC <- data[grepl(",", data[,"Chromosome"]),]
data <- data[!grepl(",", data[,"Chromosome"]),]

sepRows <- function(x)
{
	tmpC <- x["Chromosome"][[1]]
	tmpCA <- strsplit(tmpC, split=",")[[1]]

	tmpXAll <- x;
	tmpXAll["Chromosome"] <- tmpCA[1];
	tmpXAll["Target_Position"] <- tmpCA[1];

	for(i in 2:length(tmpCA))
	{
		tmpX <- x;
		tmpX["Chromosome"] <- tmpCA[i];
		tmpX["Target_Position"] <- tmpCA[i];
		tmpXAll <- rbind(tmpXAll, tmpX)
	}
	return(tmpXAll);
}
dataMC <- do.call("rbind", apply(dataMC, FUN=sepRows, MARGIN=1));
rownames(dataMC) <- c(1:nrow(dataMC))
dataMC <- data.frame(dataMC);
data <- rbind(data, dataMC); # Now 5083 rows

# 2) Now let's remove any duplicate rows
data[,2] <- factor(data[,2], levels=c("BrainTumorV4_LFS", "NonCNS_Solid_v4_LFS", "Heme_nonmalignant_SPM", "LiquidTumors_v5_SPM", "braintumor", "hemenm_spm", "noncns_solid", "liquidtumor"))
data <- data[order(data[,2]),]
tmpData <- data[,c("Sample_ID", "X3_Transcript", "X5_Transcript", "Amino_Acid_Change", "CNV_Type", "Chromosome_Band", "Chromosome", "Type_varType", "Variant_Tier", "Gene")]
data <- data[!duplicated(tmpData),] # 4977 rows 

data <- data[!data[,"Variant_Tier"]=="do not report",]
write.table(data, "../data/CleanDataFinal_V1.txt", sep="\t", row.names=F)
















