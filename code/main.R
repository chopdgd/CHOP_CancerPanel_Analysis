########################################
#Author: Pichai Raman, Komal Rathi
#Purpose: Summarize data
#Date: 02/01/2018
########################################

# Call libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(randomcoloR)
source("pubTheme.R")
source("proteinCodes.R")

#Read in cleaned data
data <- read.delim("../data/CleanDataFinal_V9.txt")

########################################
#List of all Figures for paper
#Plot 1. Bar Chart of tumor types & samples
#Plot 2a. Oncoprint of mutations and lesions across Brain
#Plot 2b. Oncoprint of mutations and lesions across Heme 
#Plot 2c. Oncoprint of mutations and lesions across Liquid
#Plot 2d. Oncoprint of mutations and lesions across Solid
#Plot 3. Bar chart of mutations by Variant
#Plot 4. Stacked bar chart of Cases by Tier
########################################

# Plot 1A
# looks like 413 cases in all and 315 of just the 
tmpDat <- unique(data[,c("Sample_ID", "Type")])
tmpDat <- table(tmpDat[,"Type"])
tmpDat <- data.frame(t(tmpDat))[2:3]
p <- ggplot(tmpDat, aes(Var2, Freq))+geom_bar(stat="identity")+theme_Publication()
p <- p+xlab("Cancer Type")+ylab("Frequency")+ggtitle(paste("Sample Frequency (n=", sum(tmpDat[,2]), ")", sep=""))
p
ggsave("../figures/08222018/CancerTypeFrequencyBarChart.eps", width=9, height=5)

# remove Tier 3 fusion from high tiers
highTiers <- c("1A", "1B", "2","Tier 1B", "CNV Tier 1", "CNV Tier 1A", "CNV Tier 1B", "CNV Tier 2", "Fusion Tier 1", "Fusion Tier 2", "Fusion Tier 1A")
data$Variant_Tier <- as.character(data$Variant_Tier)
dataHT <- data[data$Variant_Tier %in% highTiers,]

# Plot 1B - only cases with at least one of these mutations
tmpDat <- unique(dataHT[,c("Sample_ID", "Type")])
tmpDat <- table(tmpDat[,"Type"])
tmpDat <- data.frame(t(tmpDat))[2:3]
p <- ggplot(tmpDat, aes(Var2, Freq))+geom_bar(stat="identity")+theme_Publication()
p <- p+xlab("Cancer Type")+ylab("Frequency")+ggtitle(paste("Sample Frequency (n=", sum(tmpDat[,2]), ")", sep=""))
p
ggsave("../figures/08222018/CancerTypeFrequencyBarChartHighTier.eps", width=9, height=5)

# Plot 2, 3, 4
# Function
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  }, 
  Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "red", col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)
col = c(Mutation = "red", Fusion="#008000")

createPlot2 <- function(dataHT, x, freqCut=1) { 
  #Create Oncoprint
  
  #Pull out rows corresponding to panel type
  tmpDatOrig <- dataHT[dataHT[,"Type"]==x,]
  tmpDatOrig[,"Gene"] <- as.character(tmpDatOrig[,"Gene"])
  tmpDatOrig[,"Sample_ID"] <- as.character(tmpDatOrig[,"Sample_ID"])
  
  #Remove CNVs
  tmpDatOrig <- tmpDatOrig[!grepl("CNV", tmpDatOrig[,"Variant_Tier"]), ]
  
  #Get to pertinent columns cha is just for debugging
  tmpDat <- unique(tmpDatOrig[,c("Sample_ID", "Gene")])
  cha <- unique(tmpDatOrig[,c("Sample_ID", "Gene", "Subcategory")])
  
  #Filter genes that aren't mutated frequently
  filtGenes <- names(which(table(tmpDat[,"Gene"])>freqCut))
  tmpDat <- tmpDat[tmpDat[,"Gene"]%in%filtGenes,]
  
  tmpDat[,"Value"] <- ifelse(grepl("-", tmpDat[,"Gene"]), "Fusion", "Mutation")
  colnames(tmpDat)[1] <- "sampleID"
  tmpDat[,1] <- as.character(tmpDat[,1])
  tmpDat <- tmpDat %>% spread(Gene, Value, fill="")
  rownames(tmpDat) <- tmpDat[,1]
  tmpDat <- tmpDat[-1]
  tmpDat <- t(tmpDat)
  
  #Create Bottom Annotation
  tmpDatAnnot <- unique(tmpDatOrig[,c("Sample_ID", "Subcategory", "Age_Years")])
  tmpDatAnnot[,"Subcategory"]<- as.character(tmpDatAnnot[,"Subcategory"])
  rownames(tmpDatAnnot) <- tmpDatAnnot[,1]
  tmpDatAnnot <- tmpDatAnnot[colnames(tmpDat),]
  tmpDatAnnot <- tmpDatAnnot[-1]
  subCat <- unique(as.character(tmpDatAnnot[,"Subcategory"]))
  subCat <- subCat[!is.na(subCat)]
  cols <- c("#A6CEE3", "#1F78B4", "#000080", "#B2DF8A", "#33A02C", "#004E00","#FF8DA1","#E31A1C","#FFEA00","#FDBF6F","#FF7F00","#7C3720","#CAB2D6","#6A3D9A")
  subCatCols <- cols[1:length(subCat)]
  names(subCatCols)<- subCat
  ageYears <- tmpDatAnnot[,"Age_Years"]
  ha1 = HeatmapAnnotation(df=tmpDatAnnot, 
                          col = list(Subcategory = subCatCols,
                                     Age_Years = colorRamp2(c(min(ageYears, na.rm=T), max(ageYears, na.rm=T)), c("white", "red"))),
                          which="column")
  oncoPrint(tmpDat, alter_fun=alter_fun, col=col,
            remove_empty_columns = F,
            column_title = paste(x, "Oncoprint", sep=" "),
            bottom_annotation = ha1
  )
}

# remove NUDT51, TPMT, MPL from Brain
png("../figures/08222018/OncoprintBrainTumor.png", height=1440, width=2880, res=200)
createPlot2(dataHT = dataHT[-which(dataHT$Gene %in% c("NUDT15","TPMT","MPL")),], "Brain Tumor")
dev.off()

# keep MPL here
png("../figures/08222018/OncoprintLiquidTumors.png", height=1440, width=2880, res=200)
createPlot2(dataHT = dataHT[-which(dataHT$Gene %in% c("NUDT15","TPMT")),], "Liquid Tumor")
dev.off()

# remove all three from Solid tumors
png("../figures/08222018/OncoprintNonCNS_Solid.png", height=1440, width=2880, res=200)
createPlot2(dataHT = dataHT[-which(dataHT$Gene %in% c("NUDT15","TPMT","MPL")),], "Solid Tumor (Non-CNS)")
dev.off()

# Plot 5 
tmpData <- dataHT[-which(dataHT$Gene %in% c("NUDT15","TPMT")),c("Type", "Gene", "Protein_Change")]
tmpData <- tmpData[!tmpData[,3]=="",]
for (i in 1:length(code3)){
  tmpData$Protein_Change <- gsub(code3[i], code1[i], tmpData$Protein_Change, ignore.case=TRUE)
}
tmpData[,"Mutation"] <- paste(tmpData[,"Gene"], tmpData[,"Protein_Change"])
tmpData <- tmpData[,c("Type", "Mutation")]
tmpData <- table(tmpData)
tmpData <- as.data.frame(as.matrix(tmpData))
tmpData <- tmpData[tmpData$Freq>1,]
ct <- plyr::count(tmpData, 'Mutation', 'Freq')
tmpData <- merge(tmpData, ct, by = 'Mutation')
tmpData <- tmpData[which(tmpData$Mutation != "MPL p.K39N"),]
tmpData$Mutation <- reorder(tmpData$Mutation, tmpData$freq)
p <- ggplot(tmpData, aes(Mutation, Freq, fill = Type)) + 
  geom_bar(stat="identity") + 
  theme_Publication() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)) + 
  scale_y_continuous(breaks= pretty_breaks(n = 10))
png("../figures/08222018/ProteinChange.png", height=1440, width=2880, res=200)
p
dev.off()

# Plot5 redo: merge x-axis to protein level
tmpData <- dataHT[-which(dataHT$Gene %in% c("NUDT15","TPMT")),c("Type", "Gene", "Protein_Change")]
tmpData <- tmpData[!tmpData[,3]=="",]
for (i in 1:length(code3)){
  tmpData$Protein_Change <- gsub(code3[i], code1[i], tmpData$Protein_Change, ignore.case=TRUE)
}
tmpData <- tmpData[which(tmpData$Protein_Change != "p.K39N"),]
tmpData[,"Mutation"] <- paste(tmpData[,"Gene"], tmpData[,"Protein_Change"])
tmpData <- tmpData[,c("Type", "Mutation")]
tmpData <- table(tmpData)
tmpData <- as.data.frame(as.matrix(tmpData))
tmpData <- cbind(tmpData, colsplit(tmpData$Mutation, pattern = ' ', names = c('Genes','Protein.change')))
ct <- plyr::count(tmpData, 'Genes', 'Freq')
tmpData <- merge(tmpData, ct, by = 'Genes')
#tmpData <- tmpData[which(tmpData$Genes != "MPL"),]
tmpData$Genes <- reorder(tmpData$Genes, tmpData$freq)
p <- ggplot(tmpData, aes(Genes, Freq, fill = Type)) + 
  geom_bar(stat="identity") + 
  theme_Publication() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)) + 
  scale_y_continuous(breaks= pretty_breaks(n = 10))
png("../figures/08222018/ProteinChange_genelevel.png", height=1240, width=3700, res=200)
p
dev.off()

# Supplementary Figure 1
lowTiers <- c("3", "4", "CNV Tier 3") 
dataLT <- data[data[,"Variant_Tier"] %in% lowTiers,]
tmpData <- dataLT[,c("Type", "Gene", "Protein_Change")]
tmpData <- tmpData[!tmpData[,3] == "",]
# # For each code replace 3letter code by 1letter code:
for (i in 1:length(code3)){
  tmpData$Protein_Change <- gsub(code3[i], code1[i], tmpData$Protein_Change, ignore.case=TRUE)
}
tmpData[,"Mutation"] <- paste(tmpData[,"Gene"], tmpData[,"Protein_Change"])
tmpData <- tmpData[,c("Type", "Mutation")]
tmpData <- table(tmpData)
tmpData <- as.data.frame(as.matrix(tmpData))
tmpData <- tmpData[tmpData[,"Freq"]>2,]
ct <- plyr::count(tmpData, 'Mutation', 'Freq')
tmpData <- merge(tmpData, ct, by = 'Mutation')
tmpData$Mutation <- reorder(tmpData$Mutation, tmpData$freq)
p <- ggplot(tmpData, aes(Mutation, Freq, fill = Type)) + 
  geom_bar(stat="identity") + theme_Publication() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)) + 
  scale_y_continuous(breaks= pretty_breaks())
png("../figures/08222018/ProteinChangeLowTiers.png", height=1440, width=2880, res=200)
p
dev.off()

# Plot 6 
tmpData <- unique(data[,c("Sample_ID", "Type", "Variant_Tier")])
tmpData <- tmpData[!tmpData$Variant_Tier=="",c("Type","Variant_Tier")]
tmpData <- table(tmpData)
tmpData <- as.data.frame(as.matrix(tmpData))
tmpData <- tmpData[tmpData$Freq>0,]
p <- ggplot(tmpData, aes(Variant_Tier, Freq)) + geom_bar(stat="identity") + facet_wrap(~Type, nrow = 3) + theme_Publication()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
png("../figures/08222018/FreqLesionsByVariantTier.png", height=1440, width=1000, res=200)
p
dev.off()



















