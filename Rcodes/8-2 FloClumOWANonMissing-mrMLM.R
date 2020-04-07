############# this is prepare the dataset for mrMLMM software
############# this is prepare the dataset for mrMLMM software
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
allsnp2sub <- read.csv("data/allsnp2sub.csv")

###Trait for Flowering traits
###Trait for Flowering traits
allsnp2sub <- read.csv("data/allsnp2sub.csv")
##import the SNPs files
##import the SNPs files
myGMmrMLMM <- read.csv("data/myGMf.96.4814.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDf.96.4814.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMf.96.4814.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYf.96.4814.csv")
str(myYmrMlMM)

myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
names(myYmrMlMM)
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMf.96.4814.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMf.96.4814.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMf.96.4814.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:24,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/F119")

####f.106.4185
####f.106.4185
myGMmrMLMM <- read.csv("data/myGMf.106.4185.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDf.106.4185.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMf.106.4185.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYf.106.4185.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
names(myYmrMlMM)
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
names(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMf.106.4185.csv", row.names = FALSE, na = "NA")

###import the myQ
myQ<- read.csv("data/myQf.106.4185.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/subgenomrMLMMflo.csv
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMf.106.4185.csv")
str(mymrMLMMflo)
mymrMLMMQflo <- myQ[myQ$Taxa %in% mymrMLMMflo$X.phenotype., ]
str(mymrMLMMQflo)
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "Taxa")] <- "<Trait>"
str(mymrMLMMQflo)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQflof.106.4185.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQflof.106.4185.csv",header=F)
str(mymrMLMMQflo)
##
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQflof.106.4185.1.csv",row.names = FALSE)


library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMf.106.4185.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMf.106.4185.csv",
      fileKin=NULL,filePS="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mymrMLMMQflof.106.4185.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:24,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=TRUE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/F106.1")

####f.116.3098
####f.116.3098
myGMmrMLMM <- read.csv("data/myGMf.116.3098.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDf.116.3098.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMf.116.3098.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYf.116.3098.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMf.116.3098.csv", row.names = FALSE, na = "NA")

###import the myQ
myQ<- read.csv("data/myQf.116.3098.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/subgenomrMLMMflo.csv
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMf.116.3098.csv")
str(mymrMLMMflo)
mymrMLMMQflo <- myQ[myQ$Taxa %in% mymrMLMMflo$X.phenotype., ]
str(mymrMLMMQflo)
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "Taxa")] <- "<Trait>"
str(mymrMLMMQflo)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQflof.116.3098.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQflof.116.3098.csv",header=F)
str(mymrMLMMQflo)
##
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQflof.116.3098.1.csv",row.names = FALSE)

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMf.116.3098.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMf.116.3098.csv",
      fileKin=NULL,filePS="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mymrMLMMQflof.116.3098.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:24,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/F116UP2")


###f.96.4814
###f.96.4814
myGMmrMLMM <- read.csv("data/myGMf.96.4814.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDf.96.4814.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMf.96.4814.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYf.96.4814.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMf.96.4814.csv", row.names = FALSE, na = "NA")
###import the myQ
myQ<- read.csv("data/myQf.96.4814.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/subgenomrMLMMflo.csv
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMf.96.4814.csv")
str(mymrMLMMflo)
mymrMLMMQflo <- myQ[myQ$Taxa %in% mymrMLMMflo$X.phenotype., ]
str(mymrMLMMQflo)
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "Taxa")] <- "<Trait>"
str(mymrMLMMQflo)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQflof.96.4814.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQflof.96.4814.csv",header=F)
str(mymrMLMMQflo)
##
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQflof.96.4814.1.csv",row.names = FALSE)

install.packages("mrMLM")
install.packages("ggplot2")
library(ggplo2)
library(mrMLM)
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMf.96.4814.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMf.96.4814.csv",
      fileKin=NULL,filePS="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mymrMLMMQflof.96.4814.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:24,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/F122UP2")

install.packages("mrMLM.GUI")
library(mrMLM.GUI)
mrMLM.GUI()



####C.124.2560
####C.124.2560
myGMmrMLMM <- read.csv("data/myGMC.124.2560.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.124.2560.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.124.2560.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYC.124.2560.csv")
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMC.124.2560.csv", row.names = FALSE, na = "NA")

###import the myQ
myQ<- read.csv("data/myQC.124.2560.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/subgenomrMLMMflo.csv
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMC.124.2560.csv")
str(mymrMLMMflo)
mymrMLMMQflo <- myQ[myQ$Taxa %in% mymrMLMMflo$X.phenotype., ]
str(mymrMLMMQflo)
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "Taxa")] <- "<Trait>"
str(mymrMLMMQflo)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQC.124.2560.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQC.124.2560.csv",header=F)
str(mymrMLMMQflo)
##
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQC.124.2560.1.csv",row.names = FALSE)


library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.124.2560.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.124.2560.csv",
      fileKin=NULL,filePS="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mymrMLMMQC.124.2560.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=13,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/C125.1")


####C.116.3293
####C.116.3293
myGMmrMLMM <- read.csv("data/myGMC.116.3293.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.116.3293.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.116.3293.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYC.116.3293.csv")
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMC.116.3293.csv", row.names = FALSE, na = "NA")

###import the myQ
myQ<- read.csv("data/myQC.116.3293.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/subgenomrMLMMflo.csv
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMC.116.3293.csv")
str(mymrMLMMflo)
mymrMLMMQflo <- myQ[myQ$Taxa %in% mymrMLMMflo$X.phenotype., ]
str(mymrMLMMQflo)
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "Taxa")] <- "<Trait>"
str(mymrMLMMQflo)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQC.116.3293.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQC.116.3293.csv",header=F)
str(mymrMLMMQflo)
##
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQC.116.3293.1.csv",row.names = FALSE)

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.116.3293.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.116.3293.csv",
      fileKin=NULL,filePS="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mymrMLMMQC.116.3293.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=13,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/C116.1")




####C.106.4202
####C.106.4202
myGMmrMLMM <- read.csv("data/myGMC.106.4202.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.106.4202.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.106.4202.csv", row.names = FALSE, na = "NA")
str(subgenomrMLMM)

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYC.106.4202.csv")
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
names(myYmrMlMM)
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMC.106.4202.csv", row.names = FALSE, na = "NA")

###import the myQ
myQ<- read.csv("data/myQC.106.4202.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/subgenomrMLMMflo.csv
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMC.106.4202.csv")
str(mymrMLMMflo)
mymrMLMMQflo <- myQ[myQ$Taxa %in% mymrMLMMflo$X.phenotype., ]
str(mymrMLMMQflo)
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "Taxa")] <- "<Trait>"
str(mymrMLMMQflo)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQC.106.4202.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQC.106.4202.csv",header=F)
str(mymrMLMMQflo)
##
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQC.106.4202.1.csv",row.names = FALSE)

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.106.4202.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.106.4202.csv",
      fileKin=NULL,filePS="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mymrMLMMQC.106.4202.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=13,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/C106.1")
library(mrMLM.GUI)
mrMLM.GUI()


#####OWA
#####OWA
#####OWA
###for OWA with 102.4322.csv
###for OWA with 102.4322.csv
####import all of the genotype data
myY <- read.csv("data/myYO.102.4322.csv")
myY <- myY[,c(1,3:4)]
myY$Taxa <- make.names(myY$Taxa)
colnames(myY)[which(names(myY) == "Taxa")] <- "<phenotype>"
source("Function/Changecolname.R")
myY <- colRename(myY)
str(myY)
write.csv(myY, file = "mrMLMM2/mrmyYO.102.4322.csv", row.names = FALSE, na = "NA")
myGD <- read.csv("data/myGDO.102.4322.csv")
str(myGD)
myGD <- data.frame(myGD,row.names = 1)
myGDtran <- data.frame(t(myGD))-1
str(myGDtran)
source("Function/FirstColumn.R")
myGDtran <- FirstColumn(myGDtran)
str(myGDtran)
colnames(myGDtran)[which(names(myGDtran) == "Taxa")] <- "rs"
myGM <- read.csv("data/myGMO.102.4322.csv")
str(myGM)
colnames(myGM)[which(names(myGM) == "Name")] <- "rs"
str(myGM)
myQ<- read.csv("data/myQO.102.4322.csv")
allsnp1 <- read.csv("data/allsnp1.csv")
str(allsnp1)
mrMLMMgeno <-allsnp1[match(colnames(myGD), allsnp1$rs, nomatch=0),]
mrMLMMgeno <- mrMLMMgeno[,c(2,5)]
str(mrMLMMgeno)
mrMLMMgeno2 <- plyr::join_all(list(data.frame(myGM),mrMLMMgeno,myGDtran),by="rs")
str(mrMLMMgeno2)
###change the name in order to fit the software requirment
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "rs")] <- "rs#"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Chromosome")] <- "chrom"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Position")] <- "pos"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "genotype.for.code.1")] <- "genotype for code 1"
str(mrMLMMgeno2)
write.csv(mrMLMMgeno2, file = "mrMLMM2/genomrmyYO.102.4322.csv", row.names = FALSE, na = "NA")

###import the myQ
myQ<- read.csv("data/myQO.102.4322.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/myYmrMLMMculm3o.csv
mymrMLMMflo <- read.csv("mrMLMM2/mrmyYO.102.4322.csv")
str(mymrMLMMflo)
mymrMLMMQflo <- myQ[myQ$Taxa %in% mymrMLMMflo$X.phenotype., ]
str(mymrMLMMQflo)
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "Taxa")] <- "<Trait>"
str(mymrMLMMQflo)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQO.102.4322.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQO.102.4322.csv",header=F)
str(mymrMLMMQflo)
##
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQO.102.4322.1.csv",row.names = FALSE)


install.packages("mrMLM")
install.packages("tibble")
install.packages("zoo")
install.packages("lpSolve")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\genomrmyYO.102.4322.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mrmyYO.102.4322.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/OWA106")


####import all of the genotype data
myY <- read.csv("data/myYO.112.3450.csv")
myY <- myY[,c(1,3:4)]
myY$Taxa <- make.names(myY$Taxa)
colnames(myY)[which(names(myY) == "Taxa")] <- "<phenotype>"
source("Function/Changecolname.R")
myY <- colRename(myY)
str(myY)
write.csv(myY, file = "mrMLMM2/mrmyYO.112.3450.csv", row.names = FALSE, na = "NA")
myGD <- read.csv("data/myGDO.112.3450.csv")
str(myGD)
myGD <- data.frame(myGD,row.names = 1)
myGDtran <- data.frame(t(myGD))-1
str(myGDtran)
source("Function/FirstColumn.R")
myGDtran <- FirstColumn(myGDtran)
str(myGDtran)
colnames(myGDtran)[which(names(myGDtran) == "Taxa")] <- "rs"
myGM <- read.csv("data/myGMO.112.3450.csv")
str(myGM)
colnames(myGM)[which(names(myGM) == "Name")] <- "rs"
str(myGM)
myQ<- read.csv("data/myQO.112.3450.csv")
allsnp1 <- read.csv("data/allsnp1.csv")
str(allsnp1)
mrMLMMgeno <-allsnp1[match(colnames(myGD), allsnp1$rs, nomatch=0),]
mrMLMMgeno <- mrMLMMgeno[,c(2,5)]
str(mrMLMMgeno)
mrMLMMgeno2 <- plyr::join_all(list(data.frame(myGM),mrMLMMgeno,myGDtran),by="rs")
str(mrMLMMgeno2)
###change the name in order to fit the software requirment
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "rs")] <- "rs#"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Chromosome")] <- "chrom"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Position")] <- "pos"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "genotype.for.code.1")] <- "genotype for code 1"
str(mrMLMMgeno2)
write.csv(mrMLMMgeno2, file = "mrMLMM2/genomrmyYO.112.3450.csv", row.names = FALSE, na = "NA")



###import the myQ
myQ<- read.csv("data/myQO.112.3450.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/myYmrMLMMculm3o.csv
mymrMLMMflo <- read.csv("mrMLMM2/mrmyYO.112.3450.csv")
str(mymrMLMMflo)
mymrMLMMQflo <- myQ[myQ$Taxa %in% mymrMLMMflo$X.phenotype., ]
str(mymrMLMMQflo)
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "Taxa")] <- "<Trait>"
str(mymrMLMMQflo)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQO.112.3450.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQO.112.3450.csv",header=F)
str(mymrMLMMQflo)
##
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQO.112.3450.1.csv",row.names = FALSE)

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\genomrmyYO.112.3450.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mrmyYO.112.3450.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/OWA116")


####import all of the genotype data
myY <- read.csv("data/myYO.122.2646.csv")
#myY <- myY[,c(1,3:4)]
myY$Taxa <- make.names(myY$Taxa)
colnames(myY)[which(names(myY) == "Taxa")] <- "<phenotype>"
source("Function/Changecolname.R")
myY <- colRename(myY)
str(myY)
write.csv(myY, file = "mrMLMM2/mrmyYO.122.2646.csv", row.names = FALSE, na = "NA")
myGD <- read.csv("data/myGDO.122.2646.csv")
str(myGD)
myGD <- data.frame(myGD,row.names = 1)
myGDtran <- data.frame(t(myGD))-1
str(myGDtran)
source("Function/FirstColumn.R")
myGDtran <- FirstColumn(myGDtran)
str(myGDtran)
colnames(myGDtran)[which(names(myGDtran) == "Taxa")] <- "rs"
myGM <- read.csv("data/myGMO.122.2646.csv")
str(myGM)
colnames(myGM)[which(names(myGM) == "Name")] <- "rs"
str(myGM)
myQ<- read.csv("data/myQO.122.2646.csv")
allsnp1 <- read.csv("data/allsnp1.csv")
str(allsnp1)
mrMLMMgeno <-allsnp1[match(colnames(myGD), allsnp1$rs, nomatch=0),]
mrMLMMgeno <- mrMLMMgeno[,c(2,5)]
str(mrMLMMgeno)
mrMLMMgeno2 <- plyr::join_all(list(data.frame(myGM),mrMLMMgeno,myGDtran),by="rs")
str(mrMLMMgeno2)
###change the name in order to fit the software requirment
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "rs")] <- "rs#"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Chromosome")] <- "chrom"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Position")] <- "pos"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "genotype.for.code.1")] <- "genotype for code 1"
str(mrMLMMgeno2)
write.csv(mrMLMMgeno2, file = "mrMLMM2/genomrmyYO.122.2646.csv", row.names = FALSE, na = "NA")
###import the myQ
myQ<- read.csv("data/myQO.122.2646.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/myYmrMLMMculm3o.csv
mymrMLMMflo <- read.csv("mrMLMM2/mrmyYO.122.2646.csv")
str(mymrMLMMflo)
mymrMLMMQflo <- myQ[myQ$Taxa %in% mymrMLMMflo$X.phenotype., ]
str(mymrMLMMQflo)
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "Taxa")] <- "<Trait>"
str(mymrMLMMQflo)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
#write.csv(mymrMLMMQflo,file ="data/mymrMLMMQO.122.2646.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQO.122.2646.csv",header=F)
str(mymrMLMMQflo)
##
colnames(mymrMLMMQflo)[which(names(mymrMLMMQflo) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQflo <- colRename(mymrMLMMQflo)
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQO.122.2646.1.csv",row.names = FALSE)


install.packages("mrMLM")
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\genomrmyYO.122.2646.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mrmyYO.122.2646.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=2,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=TRUE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/O122")

####with example has PS and without PS
###example with PS
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.106.4202.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.106.4202.csv",
      fileKin=NULL,filePS="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mymrMLMMQC.106.4202.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=13,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/C106.1")
###example withou PS
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\genomrmyYO.122.2646.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mrmyYO.122.2646.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=2,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=TRUE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/O122")