####rrBLUP method
####import the genotype data with NON missing values SNPs 
###flowering f96
####import the genotype data with NON missing values SNPs 
myY <- read.csv("data/myYf.96.4814.csv")
myGD <- read.csv("data/myGDf.96.4814.csv",row.names=1)
myGM <- read.csv("data/myGMf.96.4814.csv")
myQ<- read.csv("data/myQf.96.4814.csv",row.names=1)
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
gwasResults.flo96 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                           min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.flo96, file = "rrBLUPF96UP/rrBLUP_GWAS_results.flo96.csv")
str(gwasResults.flo96)
###this one for qq plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPF96UP",
                      paste("flo96 QQ plot of VAR_", names(gwasResults.flo96[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.flo96[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPF96UP",
                      paste("flo96 Manhattan plots of _", names(gwasResults.flo96[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.flo96$Name,
                       CHR = gwasResults.flo96$Chromosome,
                       BP = gwasResults.flo96$Position,
                       P = 10 ^ -gwasResults.flo96[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

gwasResults.flo96 <- read.csv("rrBLUPF96UP/rrBLUP_GWAS_results.flo96.csv",row.names=1)
str(gwasResults.flo96)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.flo96 <- adj_P_function(gwasResults.flo96, 4, 26)
str(gwasResults.flo96)
write.csv(gwasResults.flo96, file = "rrBLUPF96UP/rrBLUP_GWAS_results.flo96.csv")



###f.116.3098
###f.116.3098
####import the genotype data with NON missing values SNPs 
myY <- read.csv("data/myYf.116.3098.csv")
myGD <- read.csv("data/myGDf.116.3098.csv", row.names=1)
myGM <- read.csv("data/myGMf.116.3098.csv")
myQ<- read.csv("data/myQf.116.3098.csv", row.names=1)
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
gwasResults.flo116 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                           min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.flo116, file = "rrBLUPF116UP/rrBLUP_GWAS_results.flo116.csv")
str(gwasResults.flo116)

###this one for qq plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPF116UP",
                      paste("flo116 QQ plot of VAR_", names(gwasResults.flo116[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.flo116[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPF116UP",
                      paste("flo116 Manhattan plots of _", names(gwasResults.flo116[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.flo116$Name,
                       CHR = gwasResults.flo116$Chromosome,
                       BP = gwasResults.flo116$Position,
                       P = 10 ^ -gwasResults.flo116[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

gwasResults.flo116 <- read.csv("rrBLUPF116UP/rrBLUP_GWAS_results.flo116.csv",row.names=1)
str(gwasResults.flo116)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.flo116 <- adj_P_function(gwasResults.flo116, 4, 26)
str(gwasResults.flo116)
write.csv(gwasResults.flo116, file = "rrBLUPF116UP/rrBLUP_GWAS_results.flo116.csv")
str(gwasResults.flo116)

