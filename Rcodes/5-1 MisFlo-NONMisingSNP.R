############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/alltraits.csv", na.strings = c("",".","NA"),row.names = 1)
#allblup <- allblup[!(allblup$Entry==7|allblup$Entry==11|allblup$Entry==129|allblup$Entry==148|allblup$Entry==149|allblup$Entry==150),]
str(allblup)
allblupflo <- allblup[,c(1, 14:16,20:39)]
str(allblupflo)
allblupflo<- allblupflo[-which(rowSums(is.na(allblupflo)) == 23),]
str(allblupflo)

allblupflo <- droplevels(allblupflo)
str(allblupflo)
allblupflo$Taxa
allblupflo <- allblupflo[!(allblupflo$Taxa=="PMS-071" | allblupflo$Taxa=="PMS-076" | allblupflo$Taxa=="PMS-074" | allblupflo$Taxa=="PMS-075"
                        | allblupflo$Taxa=="UI10-00107-Illinois"| allblupflo$Taxa=="UI10-00008-Hortico"),]

allblupflo <- droplevels(allblupflo)
allblupflo$Taxa
############# step2 
############# step2
############download the SNP data that Lindsay gave to me ################### 
############loading the genotype SNP data with missing value called datacomb3####
load("~/Documents/MiacnathusSNPinformation /160322filteredSNPs.RData")
load("C:/Users/Admin/Desktop/New folder/miscanthus study-1/Misthcanthus CCA data analysis/160322filteredSNPs.RData")
########## inspect the genotype data
#str(datacomb3) # 594 genotypes, 46,177 markers
#dimnames(datacomb3)[[1]] # this is how the genotypes are named in the SNP data

#############changing the rowname to the first columne###########################
source("Function/subsetNONmissingSNP.GAPIT.R")
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblupflo)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 38, 
### Remove the row with 38% missing value (more than 38% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.42),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]
SNPnmcoltr <- data.frame(t(SNPnmcol))
str(SNPnmcoltr)

SNP <- SNPnmcol
  SNPmeans <- colMeans(SNP, na.rm = TRUE)
  #range(SNPmeans)
  #levels(as.factor(SNPmeans))
  # determine for which SNPs 0 is the major allele
  ZeroIsMajor <- SNPmeans < 1
  # set up an empty vector to indicate which SNPs to keep
  SNPsToKeep <- logical(length(ZeroIsMajor))
  # define a minimum number of accessions with minor allele
  minWithMinor <- 3
  # figure out which to keep based on how many individuals have minor allele
  SNPsToKeep[ZeroIsMajor] <- colSums(SNP[, ZeroIsMajor] == 1 | SNP[, ZeroIsMajor] == 2, na.rm = TRUE) >= minWithMinor
  SNPsToKeep[!ZeroIsMajor] <- colSums(SNP[, !ZeroIsMajor] == 1 | SNP[, !ZeroIsMajor] == 0, na.rm = TRUE) >= minWithMinor
  # this is how I would filter based on MAF
  MAF <- SNPmeans/2
  MAF[!ZeroIsMajor] <- 1 - MAF[!ZeroIsMajor] # flip the frequencies if allele 2 is the common one
  hist(MAF) # most near zero, max 0.5
  
  SNPsToKeep1 <- SNPsToKeep
  SNP1 <- SNP[, SNPsToKeep1]
  SNP1tr <- data.frame(t(SNP1))
  different1 <- subset(SNPnmcoltr, !rownames(SNPnmcoltr) %in% row.names(SNP1tr))
  different1tr <- t(different1)
  table(different1)
  
  
  write.csv(different1,file="data/different1.csv")
  SNPsToKeep2 <- MAF > 0.01 # update SNPsToKeep; 3 ind. with minor allele AND MAF > 0.01. 35,279 SNPs.
  # subset genotype matrix 
  ###using this mymat as new genotype for the next step. 
  SNP2 <- SNP[, SNPsToKeep2]
  write.csv(SNP2, file="data/SNP2.csv")
  SNP2tr <- data.frame(t(SNP2))
  
  different2 <- subset(SNP1tr, !rownames(SNP2tr) %in% row.names(SNP1tr))
  write.csv(different2,file="data/different2.csv")
  ###transform the data back again with SNP in the row and individual in the column





############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.116.3098 <- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
myGM$Name
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.116.3098))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.116.3098 <- getmyGD(myGM,myGD.116.3098)
myGM.116.3098 <- getmyGM(myGM,myGD.116.3098)
###check if they are shared the same name in the same order
identical(as.character(myGM.116.3098$Name), colnames(myGD.116.3098))
write.csv(myGM.116.3098, file = "data/myGMf.116.3098.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.116.3098 <- FirstColumn(myGD.116.3098)
write.csv(myGD.116.3098, file = "data/myGDf.116.3098.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.116.3098 <- allblupflo[match(myGD.116.3098$Taxa,allblupflo$Taxa, nomatch=0),]
myY.116.3098 <- droplevels(myY.116.3098)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
#### write out the phenotype with correct TAXA
write.csv(myY.116.3098, file = "data/myYf.116.3098.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.116.3098$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.116.3098 <- myQ[match(myY.116.3098$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.116.3098, file = "data/myQf.116.3098.csv",row.names = FALSE)

############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/alltraits.csv", na.strings = c("",".","NA"),row.names = 1)
#allblup <- allblup[!(allblup$Entry==7|allblup$Entry==11|allblup$Entry==129|allblup$Entry==148|allblup$Entry==149|allblup$Entry==150),]
str(allblup)
allblupflo <- allblup[,c(1, 14:16,20:39)]
str(allblupflo)
allblupflo<- allblupflo[-which(rowSums(is.na(allblupflo)) == 23),]
str(allblupflo)

allblupflo <- droplevels(allblupflo)
str(allblupflo)
allblupflo$Taxa
allblupflo <- allblupflo[!(allblupflo$Taxa=="PMS-071" | allblupflo$Taxa=="PMS-076" | allblupflo$Taxa=="PMS-074" | allblupflo$Taxa=="PMS-075"
                           | allblupflo$Taxa=="UI10-00107-Illinois"| allblupflo$Taxa=="UI10-00008-Hortico"),]

allblupflo <- droplevels(allblupflo)
############# step2 
############# step2
############download the SNP data that Lindsay gave to me ################### 
############loading the genotype SNP data with missing value called datacomb3####
load("~/Documents/MiacnathusSNPinformation /160322filteredSNPs.RData")
load("C:/Users/Admin/Desktop/New folder/miscanthus study-1/Misthcanthus CCA data analysis/160322filteredSNPs.RData")
########## inspect the genotype data
#str(datacomb3) # 594 genotypes, 46,177 markers
#dimnames(datacomb3)[[1]] # this is how the genotypes are named in the SNP data

#############changing the rowname to the first columne###########################
source("Function/subsetNONmissingSNP.GAPIT.R")
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblupflo)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 31, 
### Remove the row with 31% missing value (more than 31% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.31),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.106.4185 <- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.106.4185))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.106.4185 <- getmyGD(myGM,myGD.106.4185)
myGM.106.4185 <- getmyGM(myGM,myGD.106.4185)
###check if they are shared the same name in the same order
identical(as.character(myGM.106.4185$Name), colnames(myGD.106.4185))
write.csv(myGM.106.4185, file = "data/myGMf.106.4185.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.106.4185 <- FirstColumn(myGD.106.4185)
write.csv(myGD.106.4185, file = "data/myGDf.106.4185.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.106.4185 <- allblupflo[match(myGD.106.4185$Taxa,allblupflo$Taxa, nomatch=0),]
myY.106.4185 <- droplevels(myY.106.4185)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
#### write out the phenotype with correct TAXA
write.csv(myY.106.4185, file = "data/myYf.106.4185.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.106.4185$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.106.4185 <- myQ[match(myY.106.4185$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.106.4185, file = "data/myQf.106.4185.csv",row.names = FALSE)

############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/alltraits.csv", na.strings = c("",".","NA"),row.names = 1)
#allblup <- allblup[!(allblup$Entry==7|allblup$Entry==11|allblup$Entry==129|allblup$Entry==148|allblup$Entry==149|allblup$Entry==150),]
str(allblup)
allblupflo <- allblup[,c(1, 14:16,20:39)]
str(allblupflo)
allblupflo<- allblupflo[-which(rowSums(is.na(allblupflo)) == 23),]
str(allblupflo)

allblupflo <- droplevels(allblupflo)
str(allblupflo)
allblupflo$Taxa
allblupflo <- allblupflo[!(allblupflo$Taxa=="PMS-071" | allblupflo$Taxa=="PMS-076" | allblupflo$Taxa=="PMS-074" | allblupflo$Taxa=="PMS-075"
                           | allblupflo$Taxa=="UI10-00107-Illinois"| allblupflo$Taxa=="UI10-00008-Hortico"),]

allblupflo <- droplevels(allblupflo)
############# step2 
############# step2
############download the SNP data that Lindsay gave to me ################### 
############loading the genotype SNP data with missing value called datacomb3####
load("~/Documents/MiacnathusSNPinformation /160322filteredSNPs.RData")
load("C:/Users/Admin/Desktop/New folder/miscanthus study-1/Misthcanthus CCA data analysis/160322filteredSNPs.RData")
########## inspect the genotype data
#str(datacomb3) # 594 genotypes, 46,177 markers
#dimnames(datacomb3)[[1]] # this is how the genotypes are named in the SNP data

#############changing the rowname to the first columne###########################
source("Function/subsetNONmissingSNP.GAPIT.R")
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblupflo)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 40, 
### Remove the row with 40% missing value (more than 40% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.27),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.96.4814 <- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.96.4814))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.96.4814 <- getmyGD(myGM,myGD.96.4814)
myGM.96.4814 <- getmyGM(myGM,myGD.96.4814)
###check if they are shared the same name in the same order
identical(as.character(myGM.96.4814$Name), colnames(myGD.96.4814))
write.csv(myGM.96.4814, file = "data/myGMf.96.4814.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.96.4814 <- FirstColumn(myGD.96.4814)
write.csv(myGD.96.4814, file = "data/myGDf.96.4814.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.96.4814 <- allblupflo[match(myGD.96.4814$Taxa,allblupflo$Taxa, nomatch=0),]
myY.96.4814 <- droplevels(myY.96.4814)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
#### write out the phenotype with correct TAXA
write.csv(myY.96.4814, file = "data/myYf.96.4814.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.96.4814$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.96.4814 <- myQ[match(myY.96.4814$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.96.4814, file = "data/myQf.96.4814.csv",row.names = FALSE)






