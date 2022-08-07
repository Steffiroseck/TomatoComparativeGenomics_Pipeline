# Load the necessary Packages
  library(dplyr)
  library(GO.db)

# Set the working directory

  setwd("/home/steffi/")
  data = scan('Results/Chilense.GO.IDs.final', character(), sep='\t')
  result=select(GO.db, keys=data, columns = c("TERM",'ONTOLOGY'), keytype = "GOID")
  result=data.frame(result)
  colnames(result)=c('GO', 'TERM', 'ONTOLOGY')

# Extract all GO terms with salt or salinity or water or drought
  data_results<-result %>%  filter_all(any_vars(grepl(pattern="(salt|salinity|water|drought)", .)))
  
# Remove duplicates
  dup<-data_results[!duplicated(data_results$GO), ]
  write.table(data.frame(result),"Results/chilense.GO.terms.all", quote=F, row.names=F)

# Important GO terms related to the keywords searched are written to a file
  write.table(data.frame(dup),"Results/chilense.GO.terms.salt.drought", quote=F, row.names=F, sep=",")
