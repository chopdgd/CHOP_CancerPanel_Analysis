########################################
# Author: Komal Rathi
# Purpose: Code to do a final cleanup
# Date: 01/24/18
########################################

# Call libraries
library("tidyverse")

# fix Yes and No
fixyn <- function(x){
  x <- gsub('^no', 'No', x, ignore.case = TRUE, perl = T)
  x <- gsub('^yes', 'Yes', x, ignore.case = TRUE, perl = T)
  return(x)
}

# Read clean data
data <- read.delim("../data/CleanDataFinal_V5.txt", stringsAsFactors = F) # 4,780 rows

# fix Yes, No
cols <- grep('tumor|impact|germline', colnames(data), value = T)
data[,cols] <- apply(data[,cols], MARGIN = 2, FUN = function(x) fixyn(x))
# data[,c(26:28,33:40)] <- apply(data[,c(26:28,33:40)], MARGIN = 2, FUN = function(x) fixyn(x))

# remove trailing spaces
data.cleaned <- apply(data, 2, FUN = function(x) trimws(x, which = 'both'))
data.cleaned <- as.data.frame(data.cleaned, stringsAsFactors = F)

# write out V6 cleaned up file - send to Lea
write.table(data.cleaned, file = '../data/CleanDataFinal_V6.txt', sep = "\t", row.names = F) # 4780 rows

