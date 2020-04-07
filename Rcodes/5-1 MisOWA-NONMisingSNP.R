############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
allblup <- read.csv("data/alltraitsOWA2.csv", na.strings = c("",".","NA"),row.names = 1)

str(allblup)
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
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblup)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 38, 
### Remove the row with 38% missing value (more than 38% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.33),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.112.3450 <- Select.MAF(SNPnmcol)

str(myGD.112.3450)
row.names(myGD.112.3450)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
myGM$Name
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.112.3450))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.112.3450 <- getmyGD(myGM,myGD.112.3450)
str(myGD.112.3450)
myGM.112.3450 <- getmyGM(myGM,myGD.112.3450)
myGM.112.3450$Name
str(myGM.112.3450)
###check if they are shared the same name in the same order
identical(as.character(myGM.112.3450$Name), colnames(myGD.112.3450))
write.csv(myGM.112.3450, file = "data/myGMO.112.3450.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.112.3450 <- FirstColumn(myGD.112.3450)
str(myGD.112.3450)
write.csv(myGD.112.3450, file = "data/myGDO.112.3450.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.112.3450 <- allblup[match(myGD.112.3450$Taxa,allblup$Taxa, nomatch=0),]
myY.112.3450 <- droplevels(myY.112.3450)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.112.3450)
#### write out the phenotype with correct TAXA
write.csv(myY.112.3450, file = "data/myYO.112.3450.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.112.3450$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.112.3450 <- myQ[match(myY.112.3450$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.112.3450, file = "data/myQO.112.3450.csv",row.names = FALSE)

############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
allblup <- read.csv("data/alltraitsOWA2.csv", na.strings = c("",".","NA"),row.names = 1)

str(allblup)
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
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblup)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 0.345, 
### Remove the row with 345% missing value (more than 345% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.42),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.122.2646 <- Select.MAF(SNPnmcol)


#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
myGM$Name
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.122.2646))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.122.2646 <- getmyGD(myGM,myGD.122.2646)
str(myGD.122.2646)
myGM.122.2646 <- getmyGM(myGM,myGD.122.2646)
myGM.122.2646$Name
str(myGM.122.2646)
###check if they are shared the same name in the same order
identical(as.character(myGM.122.2646$Name), colnames(myGD.122.2646))
write.csv(myGM.122.2646, file = "data/myGMO.122.2646.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.122.2646 <- FirstColumn(myGD.122.2646)
str(myGD.122.2646)
write.csv(myGD.122.2646, file = "data/myGDO.122.2646.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.122.2646 <- allblup[match(myGD.122.2646$Taxa,allblup$Taxa, nomatch=0),]
myY.122.2646 <- droplevels(myY.122.2646)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.122.2646)
#### write out the phenotype with correct TAXA
write.csv(myY.122.2646, file = "data/myYO.122.2646.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.122.2646$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.122.2646 <- myQ[match(myY.122.2646$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.122.2646, file = "data/myQO.122.2646.csv",row.names = FALSE)

############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
allblup <- read.csv("data/alltraitsOWA2.csv", na.strings = c("",".","NA"),row.names = 1)

str(allblup)
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
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblup)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 38, 
### Remove the row with 38% missing value (more than 38% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.298),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.01
source("Function/SelectMAF-GAPIT.R")
myGD.102.4322 <- Select.MAF(SNPnmcol)
#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
myGM$Name
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.102.4322))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.102.4322 <- getmyGD(myGM,myGD.102.4322)
str(myGD.102.4322)
myGM.102.4322 <- getmyGM(myGM,myGD.102.4322)
myGM.102.4322$Name
str(myGM.102.4322)
###check if they are shared the same name in the same order
identical(as.character(myGM.102.4322$Name), colnames(myGD.102.4322))
write.csv(myGM.102.4322, file = "data/myGMO.102.4322.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.102.4322 <- FirstColumn(myGD.102.4322)
str(myGD.102.4322)
write.csv(myGD.102.4322, file = "data/myGDO.102.4322.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.102.4322 <- allblup[match(myGD.102.4322$Taxa,allblup$Taxa, nomatch=0),]
myY.102.4322 <- droplevels(myY.102.4322)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.102.4322)
#### write out the phenotype with correct TAXA
write.csv(myY.102.4322, file = "data/myYO.102.4322.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.102.4322$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.102.4322 <- myQ[match(myY.102.4322$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.102.4322, file = "data/myQO.102.4322.csv",row.names = FALSE)







