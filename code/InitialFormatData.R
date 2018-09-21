########################################
#Author: Pichai Raman
#Purpose: Format Data from Mahdi
#Date: 10/02/17
########################################

#Call libraries
library("tidyverse")
dataSolid <- read.delim("../data/cancer_samples_files/noncns_solidv4_result.txt", header=F);
dataBrain <- read.delim("../data/cancer_samples_files/braintumorv4_result.txt", header=F);
dataLiquid <- read.delim("../data/cancer_samples_files/liquidtumorv4_result.txt", header=F);
dataHeme <- read.delim("../data/cancer_samples_files/hemenm_spm_result.txt", header=F);


#Not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

#Function to convert
convertShortWide <- function(x, type)
{	
	x <- x %>% filter(V3%!in%c("", "Avg_ROI_Coverage", "Total_Reads", "Number of Samples", "%Reads_Aligned", "%ROI_1x", "Adding sample to pool"))
	out <- spread(x, key="V3", value="V4", fill="")
	out[,2] <- type;
	colnames(out)[1:2] <- c("Sample_ID", "Type")
	return(out);
}
dataSolidSW <- convertShortWide(dataSolid, "noncns_solid")
dataBrainSW <- convertShortWide(dataBrain, "braintumor")
dataLiquidSW <- convertShortWide(dataLiquid, "liquidtumor")
dataHemeSW <- convertShortWide(dataHeme, "hemenm_spm")

#Get common column names and subset to those
intCols <- Reduce(intersect, list(colnames(dataSolidSW),colnames(dataBrainSW),colnames(dataLiquidSW), colnames(dataHemeSW)))
dataSolidSW <- dataSolidSW[,intCols]
dataBrainSW <- dataBrainSW[,intCols]
dataLiquidSW <- dataLiquidSW[,intCols]
dataHemeSW <- dataHemeSW[,intCols]
finalOut <- do.call(rbind, list(dataSolidSW, dataBrainSW, dataLiquidSW, dataHemeSW))
finalOut <- finalOut[,c(2,1,3:ncol(finalOut))];
finalOut <- finalOut[finalOut[,"Print On Report"]=="Y",]


write.table(finalOut, "../data/mergedDataFromRaw.txt", sep="\t", row.names=F)














