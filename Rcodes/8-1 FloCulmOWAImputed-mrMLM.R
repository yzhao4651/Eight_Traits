
###preparing all of the dataset need for mrMLM packages
###FLo with imputed SNP
###FLo with imputed SNP
###using the one with imputed number in SNPs datasets
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
str(subgenomrMLMM)
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
str(subgenomrMLMMtran)
#mymrMLMMflo <- read.csv("data/myYmrMLMMflo.csv")# this one has no any missing values 
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMflo3.csv")# this one has two more column and then has missing values 
str(mymrMLMMflo)
###get the same
genomymrMLMMflo <-subgenomrMLMMtran[match(mymrMLMMflo$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
genomymrMLMMflo <- data.frame(genomymrMLMMflo)+1
str(genomymrMLMMflo)
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMflo <- Select.MAF(genomymrMLMMflo)
genomymrMLMMflo1 <-subgenomrMLMM[1:4][match(genomymrMLMMflo$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMflo1)[which(names(genomymrMLMMflo1) == "rs.")] <- "rn"
genomymrMLMMflo <- plyr::join_all(list(genomymrMLMMflo1,genomymrMLMMflo),by="rn")
str(genomymrMLMMflo)
###change the name in order to fit the software requirment
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "rn")] <- "rs#"
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMflo, file = "mrMLMM2/subgenomrMLMMflo.csv", row.names = FALSE, na = "NA")
###import the myQ
myQ<- read.csv("data/myQimputedSNP19.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/subgenomrMLMMflo.csv
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMflo3.csv")
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
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQflo.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQflo.csv",header=F)
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
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQflo.1.csv",row.names = FALSE)

###for all Culm
###for all culm
###this one has the SNPS with imputed number not int
###this one has the SNPS with imputed number not int
####import all of the genotype data
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
###seperate two parts
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
####import all of the phenotyp data
###Clum trait
mymrMLMMculm <- read.csv("mrMLMM2/myYmrMLMMculm3.csv")
str(mymrMLMMculm)
###get the same
genomymrMLMMculm <-subgenomrMLMMtran[match(mymrMLMMculm$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMculm <- Select.MAF(genomymrMLMMculm)
##get the mateched SNP 
genomymrMLMMculm1 <-subgenomrMLMM[1:4][match(genomymrMLMMculm$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMculm1)[which(names(genomymrMLMMculm1) == "rs.")] <- "rn"
genomymrMLMMculm <- plyr::join_all(list(genomymrMLMMculm1,genomymrMLMMculm),by="rn")
str(genomymrMLMMculm)
###change the name in order to fit the software requirment
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "rn")] <- "rs#"
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMculm, file = "data/subgenomrMLMMculm.csv", row.names = FALSE, na = "NA")

###import the myQ
myQ<- read.csv("data/myQimputedSNP19.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/myYmrMLMMculm3o.csv
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMculm3.csv")
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
write.csv(mymrMLMMQflo,file ="data/mymrMLMMQculm.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("data/mymrMLMMQculm.csv",header=F)
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
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQculm.1.csv",row.names = FALSE)



###this one for OWA2
###this one for OWA2
###using the one with imputed number in SNPs datasets
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
str(subgenomrMLMM)
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
str(subgenomrMLMMtran)
#mymrMLMMflo <- read.csv("data/myYmrMLMMflo.csv")# this one has no any missing values 
mymrMLMMSuv <- read.csv("data/myYimputedSNP19OWA2.csv")
mymrMLMMSuv$Taxa <- make.names(mymrMLMMSuv$Taxa)

#mymrMLMMSuv <- read.csv("mrMLMM2/myYmrMlMMSuv.csv")
str(mymrMLMMSuv)
###get the same
genomymrMLMMSuv <-subgenomrMLMMtran[match(mymrMLMMSuv$Taxa,rownames(subgenomrMLMMtran), nomatch=0),]
genomymrMLMMSuv <- data.frame(genomymrMLMMSuv)+1
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMSuv <- Select.MAF(genomymrMLMMSuv)
genomymrMLMMSuv1 <-subgenomrMLMM[1:4][match(genomymrMLMMSuv$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMSuv1)[which(names(genomymrMLMMSuv1) == "rs.")] <- "rn"
genomymrMLMMSuv <- plyr::join_all(list(genomymrMLMMSuv1,genomymrMLMMSuv),by="rn")
str(genomymrMLMMSuv)
###change the name in order to fit the software requirment
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "rn")] <- "rs#"
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMSuv, file = "mrMLMM2/subgenomrMLMMOWA2.csv", row.names = FALSE, na = "NA")
str(genomymrMLMMSuv)

colnames(mymrMLMMSuv)[which(names(mymrMLMMSuv) == "Taxa")] <- "<phenotype>"
colnames(mymrMLMMSuv)[which(names(mymrMLMMSuv) == "OWA")] <- "Trait1"
write.csv(mymrMLMMSuv, file="mrMLMM2/myYmrMlMMOWA2.csv",row.names = FALSE,na = "NA")


###import the myQ
myQ<- read.csv("data/myQimputedSNP19OWA2.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the mrMLMM2/subgenomrMLMMflo.csv
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMlMMOWA2.csv")
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
write.csv(mymrMLMMQflo,file ="mrMLMM2/mymrMLMMQOWA2.csv",row.names = FALSE)
###import 
mymrMLMMQflo <- read.csv("mrMLMM2/mymrMLMMQOWA2.csv",header=F)
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
write.csv(mymrMLMMQflo,file="mrMLMM2/mymrMLMMQOWA2.1.csv",row.names = FALSE)

###this is the function used for run the data 
###without PS
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\data\\subgenomrMLMMflo.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=TRUE,
      Plotformat ="jpeg",Resolution="Low", dir= "C:/Users/Admin/Desktop/Miscanthus/Miscanthus/mrMLMM2/Resultfloall1.2")

### With PS
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\data\\subgenomrMLMMflo.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mymrMLMMQflo.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=TRUE,
      Plotformat ="jpeg",Resolution="Low", dir= "C:/Users/Admin/Desktop/Miscanthus/Miscanthus/mrMLMM2/Resultfloall1.1")

