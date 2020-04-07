
GetGO <- function(allSNPSGenes.trt.mtd){
  require(dplyr)
  require(tidyr)
  allSNPSGenes.trt.mtd_new <- allSNPSGenes.trt.mtd %>% tidyr::separate(GO, c("GO1","GO2","GO3","GO4","GO5","GO6","GO7","GO8","GO9"), ",")
  allSNPSGenes.trt.mtd_new.2 <- allSNPSGenes.trt.mtd_new[,-c(11:18)]
  allSNPSGenes.trt.mtd_new.2 <- droplevels(allSNPSGenes.trt.mtd_new.2)
  colnames(allSNPSGenes.trt.mtd_new.2)[colnames(allSNPSGenes.trt.mtd_new.2)=="GO1"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.2)
  allSNPSGenes.trt.mtd_new.3 <- allSNPSGenes.trt.mtd_new[,-c(10,12:18)][!is.na(allSNPSGenes.trt.mtd_new$GO2),]
  allSNPSGenes.trt.mtd_new.3 <- droplevels(allSNPSGenes.trt.mtd_new.3)
  colnames(allSNPSGenes.trt.mtd_new.3)[colnames(allSNPSGenes.trt.mtd_new.3)=="GO2"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.3)
  allSNPSGenes.trt.mtd_new.4 <- allSNPSGenes.trt.mtd_new[,-c(10:11,13:18)][!is.na(allSNPSGenes.trt.mtd_new$GO3),]
  allSNPSGenes.trt.mtd_new.4 <- droplevels(allSNPSGenes.trt.mtd_new.4)
  colnames(allSNPSGenes.trt.mtd_new.4)[colnames(allSNPSGenes.trt.mtd_new.4)=="GO3"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.4)
  allSNPSGenes.trt.mtd_new.5 <- allSNPSGenes.trt.mtd_new[,-c(10:12,14:18)][!is.na(allSNPSGenes.trt.mtd_new$GO4),]
  allSNPSGenes.trt.mtd_new.5 <- droplevels(allSNPSGenes.trt.mtd_new.5)
  colnames(allSNPSGenes.trt.mtd_new.5)[colnames(allSNPSGenes.trt.mtd_new.5)=="GO4"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.5)
  allSNPSGenes.trt.mtd_new.6 <- allSNPSGenes.trt.mtd_new[,-c(10:13,15:18)][!is.na(allSNPSGenes.trt.mtd_new$GO5),]
  allSNPSGenes.trt.mtd_new.6 <- droplevels(allSNPSGenes.trt.mtd_new.6)
  colnames(allSNPSGenes.trt.mtd_new.6)[colnames(allSNPSGenes.trt.mtd_new.6)=="GO5"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.6)
  allSNPSGenes.trt.mtd_new.7 <- allSNPSGenes.trt.mtd_new[,-c(10:14,16:18)][!is.na(allSNPSGenes.trt.mtd_new$GO6),]
  allSNPSGenes.trt.mtd_new.7 <- droplevels(allSNPSGenes.trt.mtd_new.7)
  colnames(allSNPSGenes.trt.mtd_new.7)[colnames(allSNPSGenes.trt.mtd_new.7)=="GO6"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.7)
  allSNPSGenes.trt.mtd_new.8 <- allSNPSGenes.trt.mtd_new[,-c(10:15,17:18)][!is.na(allSNPSGenes.trt.mtd_new$GO7),]
  allSNPSGenes.trt.mtd_new.8 <- droplevels(allSNPSGenes.trt.mtd_new.8)
  colnames(allSNPSGenes.trt.mtd_new.8)[colnames(allSNPSGenes.trt.mtd_new.8)=="GO7"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.8)
  allSNPSGenes.trt.mtd_new.9 <- allSNPSGenes.trt.mtd_new[,-c(10:16,18)][!is.na(allSNPSGenes.trt.mtd_new$GO8),]
  allSNPSGenes.trt.mtd_new.9 <- droplevels(allSNPSGenes.trt.mtd_new.9)
  colnames(allSNPSGenes.trt.mtd_new.9)[colnames(allSNPSGenes.trt.mtd_new.9)=="GO8"] <- "GO"
  str(allSNPSGenes.trt.mtd_new.9)
  allSNPSGenes.trt.mtd_new.10 <- allSNPSGenes.trt.mtd_new[,-c(10:17)][!is.na(allSNPSGenes.trt.mtd_new$GO9),]
  allSNPSGenes.trt.mtd_new.10 <- droplevels(allSNPSGenes.trt.mtd_new.10)
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

