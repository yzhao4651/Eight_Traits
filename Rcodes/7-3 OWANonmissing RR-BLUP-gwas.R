####rrBLUP method
####import the genotype data with NON missing values SNPs 
###O.102.4322
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYO.102.4322.csv")
str(myY)
myGD <- read.csv("data/myGDO.102.4322.csv", row.names=1)
myGM <- read.csv("data/myGMO.102.4322.csv")
myQ<- read.csv("data/myQO.102.4322.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.S102 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                     min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.S102, file = "rrBLUPO102/rrBLUP_GWAS_results.O102.csv")
str(gwasResults.S102)

###this one for qq plots
out_start=4
out_end=5
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPO102",
                      paste("O102 QQ plot of VAR_", names(gwasResults.S102[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.S102[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=5
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPO102",
                      paste("O102 Manhattan plots of _", names(gwasResults.S102[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S102$Name,
                       CHR = gwasResults.S102$Chromosome,
                       BP = gwasResults.S102$Position,
                       P = 10 ^ -gwasResults.S102[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}


gwasResults.O102 <- read.csv("rrBLUPO102/rrBLUP_GWAS_results.O102.csv",row.names=1)
str(gwasResults.O102)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.O102 <- adj_P_function(gwasResults.O102, 4, 5)

write.csv(gwasResults.O102, file = "rrBLUPO102/rrBLUP_GWAS_results.O102.csv")


###myYO.112.3450.csv
myY <- read.csv("data/myYO.112.3450.csv")
str(myY)
myGD <- read.csv("data/myGDO.112.3450.csv", row.names=1)
myGM <- read.csv("data/myGMO.112.3450.csv")
myQ<- read.csv("data/myQO.112.3450.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05, max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.S112 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                         min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.S112, file = "rrBLUPO112/rrBLUP_GWAS_results.O112.csv")
str(gwasResults.S112)

###this one for qq plots
out_start=4
out_end=5
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPO112",
                      paste("O112 QQ plot of VAR_", names(gwasResults.S112[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.S112[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=5
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPO112",
                      paste("O112 Manhattan plots of _", names(gwasResults.S112[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S112$Name,
                       CHR = gwasResults.S112$Chromosome,
                       BP = gwasResults.S112$Position,
                       P = 10 ^ -gwasResults.S112[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

gwasResults.O112 <- read.csv("rrBLUPO112/rrBLUP_GWAS_results.O112.csv",row.names=1)
str(gwasResults.O112)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.O112 <- adj_P_function(gwasResults.O112, 4, 5)

write.csv(gwasResults.O112, file = "rrBLUPO112/rrBLUP_GWAS_results.O112.csv")


###O.122.2646
myY <- read.csv("data/myYO.122.2646.csv")
str(myY)
myGD <- read.csv("data/myGDO.122.2646.csv", row.names=1)
myGM <- read.csv("data/myGMO.122.2646.csv")
myQ<- read.csv("data/myQO.122.2646.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.S122 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                         min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.S122, file = "rrBLUPO122/rrBLUP_GWAS_results.O122.csv")
str(gwasResults.S122)

###this one for qq plots
out_start=4
out_end=5
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPO122",
                      paste("O122 QQ plot of VAR_", names(gwasResults.S122[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.S122[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=5
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPO122",
                      paste("O122 Manhattan plots of _", names(gwasResults.S122[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S122$Name,
                       CHR = gwasResults.S122$Chromosome,
                       BP = gwasResults.S122$Position,
                       P = 10 ^ -gwasResults.S122[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

gwasResults.O122 <- read.csv("rrBLUPO122/rrBLUP_GWAS_results.O122.csv",row.names=1)
str(gwasResults.O122)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.O122 <- adj_P_function(gwasResults.O122, 4, 5)

write.csv(gwasResults.O122, file = "rrBLUPO122/rrBLUP_GWAS_results.O122.csv")



