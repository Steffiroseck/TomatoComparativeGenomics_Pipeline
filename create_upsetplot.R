# Load the necessary packages

library(dplyr)
library(UpSetR)

# set the current working directory
setwd("/home/steffi/")

# Read the data table with all GO.ids
data<-read.csv("Results/GO.IDs.all.species",sep="\t")
myGeneSets <- list(
 lycopersicum=data$chilense_lycopersicum,
 lycopersicoides=data$chilense_lycopersicoides,
 pennellii=data$chilense_pennellii,
 pimpinellifolium=data$chilense_pimpinellifolium,
 sitiens=data$chilense_sitiens)

# Save the plot as pdf file
pdf("Results/Upsetplot_aminoacids.pdf")
upset(fromList(myGeneSets),order.by = "freq",
     keep.order = TRUE, text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), point.size = 3.5)
dev.off()

