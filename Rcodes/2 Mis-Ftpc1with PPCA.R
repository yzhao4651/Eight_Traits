
#### PCA analysis with ppca for flowering trait
#### PCA analysis with ppca for flowering trait
##install.packages("BiocManager")
###source("https://bioconductor.org/biocLite.R") ###intall.packages("pcaMethods")
###biocLite("pcaMethods")
###library(pcaMethods)
####import the data

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")

###this one is for all of the trait
ranefvalueall<- read.csv("data/traits1718normalited5.csv",na.strings = c("",".","NA"))
#ranefvalueall<- read.csv("data/traits1718normalited4.csv",na.strings = c("",".","NA"))
ranefvalueall <- ranefvalueall[3:39]
###check data format
str(ranefvalueall)
###PC1 for flowering labels as days
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("pcaMethods") ###intall.packages("pcaMethods")
#citation("pcaMethods")
#ranefvalueall<- read.csv("~/Documents/whole traits/ranefvalueall2.csv",na.strings = c("",".","NA"))
#ranefvalueall<- read.csv("data/ranefvalueall2.csv", na.strings = c("",".","NA"))
###rename of the column name 
#colnames(ranefvalueall)[colnames(ranefvalueall)=="X"] <- "Entry"
###check data format

ranefvaluef<- ranefvalueall[!(is.na(ranefvalueall[5])) | !(is.na(ranefvalueall[6])) | !(is.na(ranefvalueall[7]))|!(is.na(ranefvalueall[8])),]
str(ranefvaluef)

library(pcaMethods)
pc <- pca(ranefvaluef[5:8], nPcs=1, method="ppca",center = TRUE)
str(pc)

fblupimputed <- data.frame(completeObs(pc))
str(fblupimputed)
flocompletepc <- cbind(ranefvaluef[1],fblupimputed)
###PC 
fprind_comp <- prcomp(flocompletepc[2:5])
###PC scree plot
pdf(paste("PC scree plot of PC1.D", 1 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(fprind_comp)
plot(fprind_comp,type="line", main=paste0("PC scree plot of PC1.D"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###get pc1 for each individual
fprind <- cbind(ranefvaluef[1:3], fprind_comp$x[,1])
names(fprind)
colnames(fprind)[colnames(fprind)=="fprind_comp$x[, 1]"] <- "fprind"


###PC1 for flowering labels as week
###PC1 for flowering labels as week
ranefvaluef<- ranefvalueall[!(is.na(ranefvalueall[26]))|!(is.na(ranefvalueall[27])) | !(is.na(ranefvalueall[28]))|!(is.na(ranefvalueall[29])),]
str(ranefvaluef)
library(pcaMethods)
pc <- pca(ranefvaluef[26:29], nPcs=1, method="ppca",center = TRUE)
fblupimputed <- data.frame(completeObs(pc))
flocompletepc <- cbind(ranefvaluef[1],fblupimputed)
###PC 
fprinW_comp <- prcomp(flocompletepc[2:5])
###PC scree plot
pdf(paste("PC scree plot of PC1.W", 2 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(fprinW_comp)
plot(fprinW_comp,type="line", main=paste0("PC scree plot of PC1.W"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###get pc1 for each individual
fprinW <- cbind(ranefvaluef[1:3], fprinW_comp$x[,1])
names(fprinW)
colnames(fprinW)[colnames(fprinW)=="fprinW_comp$x[, 1]"] <- "fprinW"

###PC1 for flowering labels as two weeks 
###PC1 for flowering labels as two weeks 
#citation("pcaMethods")
ranefvaluef<- ranefvalueall[!(is.na(ranefvalueall[30]))|!(is.na(ranefvalueall[31])) | !(is.na(ranefvalueall[32]))|!(is.na(ranefvalueall[33])),]
library(pcaMethods)
pc <- pca(ranefvaluef[30:33], nPcs=1, method="ppca",center = TRUE)

fblupimputed <- data.frame(completeObs(pc))
flocompletepc <- cbind(ranefvaluef[1],fblupimputed)
###PC 
fprinGW_comp <- prcomp(flocompletepc[2:5])
###PC scree plot
pdf(paste("PC scree plot of PC1.GW", 3 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(fprinGW_comp)
plot(fprinGW_comp,type="line", main=paste0("PC scree plot of PC1.GW"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###get pc1 for each individual
fprinGW <- cbind(ranefvaluef[1:3], fprinGW_comp$x[,1])
names(fprinGW)
colnames(fprinGW)[colnames(fprinGW)=="fprinGW_comp$x[, 1]"] <- "fprinGW"

###PC1 for flowering labels as month 
###PC1 for flowering labels as month  
#citation("pcaMethods")
ranefvaluef<- ranefvalueall[!(is.na(ranefvalueall[34]))|!(is.na(ranefvalueall[35])) | !(is.na(ranefvalueall[36]))|!(is.na(ranefvalueall[37])),]
pc <- pca(ranefvaluef[34:37], nPcs=1, method="ppca",center = TRUE)

fblupimputed <- data.frame(completeObs(pc))
flocompletepc <- cbind(ranefvaluef[1],fblupimputed)
###PC 
fprinM_comp <- prcomp(flocompletepc[2:5])
###PC scree plot
pdf(paste("PC scree plot of PC1.M", 4 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(fprinM_comp)
plot(fprinM_comp,type="line", main=paste0("PC scree plot of PC1.M"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###get pc1 for each individual
fprinM <- cbind(ranefvaluef[1:3], fprinM_comp$x[,1])
names(fprinM)
colnames(fprinM)[colnames(fprinM)=="fprinM_comp$x[, 1]"] <- "fprinM"

###combine all ofthe traits together
str(fprind)
str(fprinW)
str(fprinGW)
str(fprinM)
flowpc <- cbind(fprind,fprinW,fprinGW,fprinM)
flowpc <- flowpc[,c(1:4,8,12,16)]
str(flowpc)
#write.csv(alltraitsblup, file = "~/Documents/whole traits/alltraitsblup.csv",row.names = T)
#write.csv(alltraitsblup, file = "data/alltraitsblup.csv",row.names = T)
#write.csv(flowpc, file = "data/flowpc.csv",row.names = F)
write.csv(flowpc, file = "data/flowpc1.csv",row.names = F)
#write.csv(alltraitsblup, file = "data/alltraitsblup2.csv",row.names = F)
#write.csv(alltraitsblup, file = "data/alltraitsblup3.csv",row.names = F)
#str(fprin1)
#write.csv(fprin1, file = "~/Documents/whole traits/fprin1.csv",row.names = T)
#str(ranefvalueall)
#citation("pcaMethods")

