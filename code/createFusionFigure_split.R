########################################
#Author: Komal Rathi
#Purpose: Create Fusion Plot
#Date: 10/29/2018
########################################

# this is the main one
library(OmicCircos)
library(biomaRt)
library(reshape2)
library(ggplot2)
data("UCSC.hg19")

## regenerate db
dat    <- data.frame(UCSC.hg19[,c(1,6,7)], n1=rep("NA", nrow(UCSC.hg19)), n1=rep("NA", nrow(UCSC.hg19)));
seg.n  <- as.character(dat[,1]);
db     <- segAnglePo(seg.dat=dat, seg=seg.n);
chr.n  <- gsub("chr", "", seg.n);
chr.po <- data.frame(chr=seg.n, po=(dat[,2]+dat[,3])/2, chr.n=chr.n) 
cols1  <- rainbow(nrow(dat), alpha=0.5);

# fusion data (remove Fusion Tier 3)
data <- read.delim("../data/CleanDataFinal_V10.txt")
fusions <- data[grep("^Fusion", data$Variant_Tier),c('Patient_ID','Sample_ID','Type','X3_Transcript','X5_Transcript','Chromosome_Band','Gene','Type_varType','Variant_Tier')]
fusions <- fusions[!fusions$Variant_Tier == "Fusion Tier 3",]
fusions <- unique(fusions)
fusions <- cbind(fusions, colsplit(fusions$Gene,'-',names = c("GeneA","GeneB")))
genes <- unique(c(fusions$GeneA, fusions$GeneB))

# biomart 
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
attributes = listAttributes(grch37)
res <- getBM(attributes = c("chromosome_name", 
                            "start_position", 
                            "end_position", "hgnc_symbol", "strand"), 
             mart = grch37, values = genes, 
             filters = "hgnc_symbol")
res.genes <- unique(res[grep('^[0-9XYM]', res$chromosome_name),])
res.genes$strand <- ifelse(res.genes$strand == 1, "+", "-")

# merge (remove IGH fusions)
fs <- merge(fusions, res.genes, by.x = "GeneA", by.y = "hgnc_symbol", all.x = T)
colnames(fs)[12:15] <- c("chr_A","start_A","end_A","strand_A")
fs <- merge(fs, res.genes, by.x = "GeneB", by.y = "hgnc_symbol", all.x = T)
colnames(fs)[16:19] <- c("chr_B","start_B","end_B","strand_B")
fs <- fs[which(fs$GeneB != "IGH"),]

# for line width
ct <- plyr::count(fs, 'Gene')
ct <- merge(ct, fs, by.x = "Gene", by.y = "Gene")
fs <- ct

fs <- fs[,c("chr_A","start_A","GeneA","chr_B","start_B","GeneB","Sample_ID","Type","freq")] 

# filter out dups
up_gene <- unique(fs$GeneA)
dw_gene <- unique(fs$GeneB)
intersect(up_gene,dw_gene)
dups <- fusions[grep("RUNX1|TFE3|ASPSCR1", fusions$Gene),]
# write.table(dups, file = '../data/fusions.txt', quote = F, sep = "\t", row.names = F)

fs <- fs[which(fs$GeneA != "BRAF"),]
fs <- fs[which(fs$GeneA != "TFE3"),]
up_gene <- unique(fs$GeneA)
dw_gene <- unique(fs$GeneB)
intersect(up_gene,dw_gene) # only RUNX1 remains


### plots
pdf(file = '../figures/FusionPlot.pdf', width = 9, height = 9)
# brain tumor
x <- fs[which(fs$Type %in% "Brain Tumor"),]

# color for each type
cols <- 'lightblue'

# lw
labs <- rbind(x[,c('chr_A','start_A','GeneA','Sample_ID','Type')], setNames(x[,c('chr_B','start_B','GeneB','Sample_ID','Type')], c('chr_A','start_A','GeneA','Sample_ID','Type')))
cts <- plyr::count(labs,c('GeneA','Type'))
cts <- merge(cts, res.genes, by.x = "GeneA", by.y = 'hgnc_symbol')
cts$color <- 'lightblue'
cols2 <- cts$color
cts <- cts[,c('chromosome_name','start_position','GeneA','freq')]

# 5'3'
labs <- rbind(x[,c('chr_A','start_A','GeneA')], setNames(x[,c('chr_B','start_B','GeneB')], c('chr_A','start_A','GeneA')))
labs <- unique(labs)
labs$color <- ifelse(labs$GeneA %in% up_gene, 'red', 'blue')
labs$color[labs$GeneA == "RUNX1"] <- "black"
labs$chr_A <- factor(labs$chr_A, levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','15','16','17','18','19','21','22','X'))
labs.col <- labs[order(labs$chr_A, labs$start_A),'color']

# first
par(mar=c(.2,.2,.2,.2))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot(c(1,800), c(1,800), type = "n", axes = FALSE, xlab = "", ylab= "")
circos(R=250, type = "chr", cir = "hg19", print.chr.lab = T, W = 10, cex = 10)
circos(R=250, cir=db, W=0, mapping=dat, type="arc2", col=cols1, lwd=10);
circos(R=280,cir="hg19",W=20, mapping=labs, type="label",side="out", col=labs.col, cex=0.5)
circos(R=140, cir="hg19", W=20, mapping=x, type="link", col=cols, lwd=2)
circos(R=160,cir="hg19", W=100, mapping=cts, col.v = 4, type = "b", col = cols2, B = FALSE,scale=TRUE, cex = 4, lwd = 2)
title("Brain Tumor", line = -2)
# legend(x = "bottomright", 
#        legend = c("5'","3'","Both"), 
#        fill = c("blue","red","black"), cex = 0.7, pt.cex = 7)

# liquid tumor
x <- fs[which(fs$Type %in% "Liquid Tumor"),]

# color for each type
cols <- 'pink'

# lw
labs <- rbind(x[,c('chr_A','start_A','GeneA','Sample_ID','Type')], setNames(x[,c('chr_B','start_B','GeneB','Sample_ID','Type')], c('chr_A','start_A','GeneA','Sample_ID','Type')))
cts <- plyr::count(labs,c('GeneA','Type'))
cts <- merge(cts, res.genes, by.x = "GeneA", by.y = 'hgnc_symbol')
cts$color <- 'pink'
cols2 <- cts$color
cts <- cts[,c('chromosome_name','start_position','GeneA','freq')]

# 5'3'
labs <- rbind(x[,c('chr_A','start_A','GeneA')], setNames(x[,c('chr_B','start_B','GeneB')], c('chr_A','start_A','GeneA')))
labs <- unique(labs)
labs$color <- ifelse(labs$GeneA %in% up_gene, 'red', 'blue')
labs$color[labs$GeneA == "RUNX1"] <- "black"
labs$chr_A <- factor(labs$chr_A, levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','15','16','17','18','19','21','22','X'))
labs.col <- labs[order(labs$chr_A, labs$start_A),'color']

# second
plot(c(1,800), c(1,800), type = "n", axes = FALSE, xlab = "", ylab= "")
circos(R=250, type = "chr", cir = "hg19", print.chr.lab = T, W = 10, cex = 10)
circos(R=250, cir=db, W=0, mapping=dat, type="arc2", col=cols1, lwd=10);
circos(R=280,cir="hg19",W=20, mapping=labs, type="label",side="out", col=labs.col, cex=0.5)
circos(R=140, cir="hg19", W=20, mapping=x, type="link", col=cols, lwd=2)
circos(R=160,cir="hg19", W=100, mapping=cts, col.v = 4, type = "b", col = cols2, B = FALSE,scale=TRUE, cex = 4, lwd = 2)
title("Liquid Tumor", line = -2)
# legend(x = "bottomright", 
#        legend = c("5'","3'","Both"), 
#        fill = c("blue","red","black"), cex = 0.7, pt.cex = 7)

# solid tumor
x <- fs[which(fs$Type %in% "Solid Tumor (Non-CNS)"),]

# color for each type
cols <- 'purple'

# lw
labs <- rbind(x[,c('chr_A','start_A','GeneA','Sample_ID','Type')], setNames(x[,c('chr_B','start_B','GeneB','Sample_ID','Type')], c('chr_A','start_A','GeneA','Sample_ID','Type')))
cts <- plyr::count(labs,c('GeneA','Type'))
cts <- merge(cts, res.genes, by.x = "GeneA", by.y = 'hgnc_symbol')
cts$color <- 'purple'
cols2 <- cts$color
cts <- cts[,c('chromosome_name','start_position','GeneA','freq')]

# 5'3'
labs <- rbind(x[,c('chr_A','start_A','GeneA')], setNames(x[,c('chr_B','start_B','GeneB')], c('chr_A','start_A','GeneA')))
labs <- unique(labs)
labs$color <- ifelse(labs$GeneA %in% up_gene, 'red', 'blue')
labs$color[labs$GeneA == "RUNX1"] <- "black"
labs$chr_A <- factor(labs$chr_A, levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','15','16','17','18','19','21','22','X'))
labs.col <- labs[order(labs$chr_A, labs$start_A),'color']

# third
plot(c(1,800), c(1,800), type = "n", axes = FALSE, xlab = "", ylab= "")
circos(R=250, type = "chr", cir = "hg19", print.chr.lab = T, W = 10, cex = 12)
circos(R=250, cir=db, W=0, mapping=dat, type="arc2", col=cols1, lwd=10);
circos(R=280,cir="hg19",W=20, mapping=labs, type="label",side="out", col=labs.col, cex=0.5)
circos(R=140, cir="hg19", W=20, mapping=x, type="link", col=cols, lwd=2)
circos(R=160,cir="hg19", W=100, mapping=cts, col.v = 4, type = "b", col = cols2, B = FALSE,scale=TRUE, cex = 4, lwd = 2)
title("Solid Tumor (Non-CNS)", line = -2)

legend(1000, 400, legend = c("5'","3'","Both"),
       fill = c("blue","red","black"),
       cex = 1, 
       title = "Subgroup",xpd = "NA")
legend(1000, 600,
       legend = c("Liquid Tumor","Brain Tumor","Solid Tumor (Non-CNS)"),
       fill = c("lightblue","pink","purple"),
       lty = c(1,1,1),
       cex = 1, 
       title = "Tumor Type", xpd = "NA")
dev.off()
