

GetGO <- function(allSNPSGenes.trt.mtd_new){
  require(dplyr)
  require(tidyr)
  allSNPSGenes.trt.mtd_new <- allSNPSGenes.trt.mtd %>% tidyr::separate(GO, c("GO1","GO2","GO3","GO4","GO5","GO6","GO7","GO8","GO9"), ",")
  str(allSNPSGenes.trt.mtd_new)
  allSNPSGenes.trt.mtd_new.2 <- allSNPSGenes.trt.mtd_new[,-c(39:46)]
  colnames(allSNPSGenes.trt.mtd_new.2)[colnames(allSNPSGenes.trt.mtd_new.2)=="GO1"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.2)
  allSNPSGenes.trt.mtd_new.3 <- allSNPSGenes.trt.mtd_new[,-c(38,40:46)][which(allSNPSGenes.trt.mtd_new$GO2 != "NA"),]
  colnames(allSNPSGenes.trt.mtd_new.3)[colnames(allSNPSGenes.trt.mtd_new.3)=="GO2"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.3)
  allSNPSGenes.trt.mtd_new.4 <- allSNPSGenes.trt.mtd_new[,-c(38:39,41:46)][which(allSNPSGenes.trt.mtd_new$GO3 != "NA"),]
  colnames(allSNPSGenes.trt.mtd_new.4)[colnames(allSNPSGenes.trt.mtd_new.4)=="GO3"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.4)
  allSNPSGenes.trt.mtd_new.5 <- allSNPSGenes.trt.mtd_new[,-c(38:40,42:46)][which(allSNPSGenes.trt.mtd_new$GO4 != "NA"),]
  colnames(allSNPSGenes.trt.mtd_new.5)[colnames(allSNPSGenes.trt.mtd_new.5)=="GO4"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.5)
  allSNPSGenes.trt.mtd_new.6 <- allSNPSGenes.trt.mtd_new[,-c(38:41,43:46)][which(allSNPSGenes.trt.mtd_new$GO5 != "NA"),]
  colnames(allSNPSGenes.trt.mtd_new.6)[colnames(allSNPSGenes.trt.mtd_new.6)=="GO5"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.6)
  allSNPSGenes.trt.mtd_new.7 <- allSNPSGenes.trt.mtd_new[,-c(38:42,44:46)][which(allSNPSGenes.trt.mtd_new$GO6 != "NA"),]
  colnames(allSNPSGenes.trt.mtd_new.7)[colnames(allSNPSGenes.trt.mtd_new.7)=="GO6"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.7)
  allSNPSGenes.trt.mtd_new.8 <- allSNPSGenes.trt.mtd_new[,-c(38:43,45:46)][which(allSNPSGenes.trt.mtd_new$GO7 != "NA"),]
  colnames(allSNPSGenes.trt.mtd_new.8)[colnames(allSNPSGenes.trt.mtd_new.8)=="GO7"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.8)
  allSNPSGenes.trt.mtd_new.9 <- allSNPSGenes.trt.mtd_new[,-c(38:44,46)][which(allSNPSGenes.trt.mtd_new$GO8 != "NA"),]
  colnames(allSNPSGenes.trt.mtd_new.9)[colnames(allSNPSGenes.trt.mtd_new.9)=="GO8"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.9)
  allSNPSGenes.trt.mtd_new.10 <- allSNPSGenes.trt.mtd_new[,-c(38:45)][which(allSNPSGenes.trt.mtd_new$GO9 != "NA"),]
  colnames(allSNPSGenes.trt.mtd_new.10)[colnames(allSNPSGenes.trt.mtd_new.10)=="GO9"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.10)
  head(allSNPSGenes.trt.mtd_new.10)
  require(lessR)
  total <- Merge(allSNPSGenes.trt.mtd_new.2, allSNPSGenes.trt.mtd_new.3)
  total <- Merge(total, allSNPSGenes.trt.mtd_new.4)
  total <- Merge(total, allSNPSGenes.trt.mtd_new.5)
  total <- Merge(total, allSNPSGenes.trt.mtd_new.6)
  total <- Merge(total, allSNPSGenes.trt.mtd_new.7)
  total <- Merge(total, allSNPSGenes.trt.mtd_new.8)
  total <- Merge(total, allSNPSGenes.trt.mtd_new.9)
  total <- Merge(total, allSNPSGenes.trt.mtd_new.10)
  return(total)
}

