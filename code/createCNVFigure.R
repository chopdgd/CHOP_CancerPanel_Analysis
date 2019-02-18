########################################
#Author: Komal Rathi
#Purpose: Create CNV Plot
#Date: 11/01/2018
########################################

################# CNV plot #####################
library(ggplot2)
library(dplyr)
source('pubTheme.R')
library(tidyr)

data <- read.delim("../data/CleanDataFinal_V10.txt")
cnv <- data[grep("CNV", data$Variant_Tier),c('Sample_ID','Type','CNV_Size','CNV_Type','Chromosome_Band','Chromosome','Type_varType','Variant_Tier')]
cnv <- unique(cnv)

# remove Y
unique(cnv$Chromosome)
cnv$Chromosome <- as.character(cnv$Chromosome)
cnv <- cnv[-which(cnv$Chromosome %in% c("X/Y","Y")),]

cnv <- cnv[which(cnv$CNV_Type %in% c("Loss","Gain")),]
cnv$CNV_x[cnv$CNV_Type == "Loss"] <- "-1"
cnv$CNV_x[cnv$CNV_Type == "Gain"] <- "1"

cnv$Chromosome_Band <- as.character(cnv$Chromosome_Band)
cnv$Chromosome_Band <- gsub('[0-9].*','',cnv$Chromosome_Band)
cnv$Chromosome_Band[cnv$Chromosome_Band == "Whole Chromosome"] <- "q"

cnv$Chromosome_Band <- ifelse(cnv$Chromosome_Band %in% c("","partial"), "p,q", cnv$Chromosome_Band)

cnv <- cnv %>% 
  mutate(Chromosome_Band = strsplit(as.character(Chromosome_Band), ",")) %>% 
  unnest(Chromosome_Band)

cnv <- cnv %>% 
  mutate(Chromosome = strsplit(as.character(Chromosome), "/")) %>% 
  unnest(Chromosome)

cnv <- unique(cnv)
cnv <- cnv[which(cnv$Chromosome != "Y"),]
cnv$label <- paste0(cnv$Chromosome,cnv$Chromosome_Band)

not.include <- c("Yp", "Yq", "13p", "14p", "15p", "18p", "20p", "21p", "22p")
cnv <- cnv[-which(cnv$label %in% not.include),]

vars <- c(1:22,"X")
vis <- c("p","q")
levs <- as.vector(t(outer(vars, vis, paste, sep="")))
cnv$label <- factor(cnv$label, levels = rev(levs))

cnv$CNV_x <- as.numeric(cnv$CNV_x)
cnv$CNV_Type <- relevel(x = cnv$CNV_Type, ref = 'Loss')

# p <- ggplot(cnv, aes(x = label, y = CNV_x, fill = Type)) + 
#   geom_bar(stat = 'identity') +
#   theme_bw() + xlab("") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_y_continuous(name="# of Patients", 
#                      labels = abs(seq(-100, 100, 10)), 
#                      breaks=seq(-100, 100, 10)) + geom_hline(yintercept = 0) + theme_Publication()
# p
# ggsave(filename = '../figures/CNV_plot.png', width = 17, height = 4, dpi = 550)

### flipped plot
# q <- ggplot(cnv, aes(x = label, y = CNV_x, fill = Type)) + 
#   geom_bar(stat = 'identity') +
#   theme_bw() + xlab("") +
#   scale_y_continuous(name="# of Patients", 
#                      labels = c('Loss', '', 'Gain'), 
#                      breaks=c(-40, 0, 40)) + geom_hline(yintercept = 0) + theme_Publication() +
#   theme(axis.ticks.x = element_blank()) + coord_flip()
# q
#ggsave(filename = '../figures/CNV_plot.png', width = 5, height = 17, dpi = 550)
################# CNV plot #####################

# first option
q <- ggplot(cnv, aes(x = label, y = CNV_x, fill = Type)) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + theme_Publication() +
  geom_hline(yintercept = 0) + 
  facet_wrap(~CNV_Type, strip.position = "bottom", scales = "free_x") +
  scale_y_continuous(expand = c(0,0),name="# of Patients", 
                     labels = c(60, 40, 20, 0, 20, 40, 60),
                     breaks=seq(-60, 60, 20)) +
  geom_text(aes(y=CNV_x * 1.1, label="")) +
  coord_flip() +
  theme(panel.spacing = unit(0, "lines")) + xlab("") + ylab("") +
  xlab("Chromosome Band") + ylab("# of Patients")
ggsave(plot = q, filename = '../figures/CNV_plot_2.png', width = 6, height = 8, dpi = 550)

# alternative
p <- ggplot(cnv, aes(x = label, y = CNV_x, fill = Type)) + 
  geom_bar(stat = 'identity', width=0.9) +
  theme_bw() + xlab("") +
  scale_y_continuous(name="# of Patients", 
                     labels = c(60, 40, 20, 0, 20, 40, 60),
                     breaks=seq(-60, 60, 20)) + 
  geom_hline(yintercept = 0, color = '#262626') + theme_Publication() +
  coord_flip()  +
  annotate("text", fontface = 4,
           size = 4,
           label = c("Loss", "Gain"), 
           x = c(20, 20), 
           y = c(-40.5, 40.5)) + xlab("Chromosome Band")
ggsave(plot = p, filename = '../figures/CNV_plot.png', width = 6, height = 8, dpi = 550)
