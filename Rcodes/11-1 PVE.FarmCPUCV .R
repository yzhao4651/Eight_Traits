### this one for MLMSUPER

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myYSRD <- read.csv("data/myYimputedSNP19SRD.csv")
myYOWA <- read.csv("data/myYimputedSNP19OWA2.csv")
myY1 <- read.csv("data/myYimputedSNP19.csv")
myY1 <- myY1[,-c(2,18)]
myY <- Reduce(function(x, y) merge(x, y, all.x=TRUE), list(myY1,myYOWA,myYSRD))

colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
levels(flo.1$Ind.SNP)
levels(flo.1$Trait.name)
flo <- flo.1[flo.1$Ind.SNP=="124+36088im" & flo.1$Method=="FarmCPU",]
str(flo)
names(flo)
flo <- droplevels(flo)
levels(flo$Ind.SNP)
levels(as.factor(flo$Method))
levels(as.factor(flo$Trait.name))
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="OWA"),]
names(flo)
names(myY)
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)

#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.15.116im")
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.25.116im")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.36088imup")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}


PVE.fun <- function(sample,myY,myGD){
  if (nrow(sample)==1){
    ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% sample$SNP])
    colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
    flowerSNPs <- merge(myY, ftdGDSNPs,by="Taxa")
    a <- ncol(myY)+1
    lm1 <-lm(flowerSNPs[[levels(sample$Trait.name)]] ~ flowerSNPs[[a]])
    af <- anova(lm1)
    afss <- af$"Sum Sq"
    PctExp <- afss/sum(afss)*100
    PVE3 <- data.frame(cbind(colnames(ftdGDSNPs[2]),PctExp[[1]]))
    colnames(PVE3) <- c("SNP","PVE")
  } else {
    ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% sample$SNP])
    colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
    flowerSNPs <- merge(myY, ftdGDSNPs,by="Taxa")
    number=1
    a <- ncol(myY)+1
    b <- a+nrow(sample)-1
    out_variable = colnames(flowerSNPs[a:b])
    outcome <- matrix(NA, nrow=ncol(flowerSNPs[a:b]),
                      ncol = 1, dimnames = list(out_variable, "PctExp"))
    for(j in a:b){
      lm1 <-lm(flowerSNPs[[levels(sample$Trait.name)]] ~ flowerSNPs[[j]])
      af <- anova(lm1)
      afss <- af$"Sum Sq"
      PctExp <- afss/sum(afss)*100
      outcome[number] <-PctExp[[1]]
      number=number+1
      PVE3 <- data.frame(outcome)
      colnames(PVE3) <-  c("PVE")
      SNP <- rownames(PVE3)
      rownames(PVE3) <- NULL
      PVE3 <- cbind(SNP,PVE3)
    }
  }
  return(PVE3)
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.36088imup/*.", extension, sep = ""))

#myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")

mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCVim <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(as.factor(FarmCPUCVim$filename))

FarmCPUCVim$filename <- gsub('FarmCPUCV0.05.36088imup/Trait.', '', FarmCPUCVim$filename)
FarmCPUCVim$filename <- gsub('.csv', '', FarmCPUCVim$filename)

levels(as.factor(FarmCPUCVim$filename))

FarmCPUCVim.1 <- unite(FarmCPUCVim,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")
flo.2$SNP.T %in% FarmCPUCVim.1$SNP.T
FarmCPUCVim.flo <- merge(flo.2,FarmCPUCVim.1,by="SNP.T")
FarmCPUCVim.flo.S <- separate(FarmCPUCVim.flo, SNP.T, c("SNP","Trait.name"),sep="/")
levels(as.factor(FarmCPUCVim.flo.S$Trait.name))
str(FarmCPUCVim.flo.S)

###for flowering 106
###for flowering 106
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYf.106.4185.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDf.106.4185.csv")
flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
flo <- flo.1[flo.1$Ind.SNP=="F106+4185" & flo.1$Method=="FarmCPU",]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- flo[!(flo$Trait.name=="Surv"),]
str(flo)
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.F106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.F106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.F106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.F106$filename)

FarmCPUCV0.05.F106$filename <- gsub('FarmCPUCV0.05.F106/Trait.', '', FarmCPUCV0.05.F106$filename)
FarmCPUCV0.05.F106$filename <- gsub('.csv', '', FarmCPUCV0.05.F106$filename)

levels(as.factor(FarmCPUCV0.05.F106$filename))

FarmCPUCV0.05.F106.1<- unite(FarmCPUCV0.05.F106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.F106.flo <- merge(flo.2,FarmCPUCV0.05.F106.1,by="SNP.T")
FarmCPUCV0.05.F106.flo.S <- separate(FarmCPUCV0.05.F106.flo, SNP.T, c("SNP","Trait.name"),sep="/")
str(FarmCPUCV0.05.F106.flo.S)

###for flowering 116
###for flowering 116
###for flowering 116
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYf.116.3098.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDf.116.3098.csv")
flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
flo <- flo.1[flo.1$Ind.SNP=="F116+3098" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.F116up")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.F116up/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.F116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.F116$filename)

FarmCPUCV0.05.F116$filename <- gsub('FarmCPUCV0.05.F116up/Trait.', '', FarmCPUCV0.05.F116$filename)
FarmCPUCV0.05.F116$filename <- gsub('.csv', '', FarmCPUCV0.05.F116$filename)

levels(as.factor(FarmCPUCV0.05.F116$filename))

FarmCPUCV0.05.F116.1 <- unite(FarmCPUCV0.05.F116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.F116.flo <- merge(flo.2,FarmCPUCV0.05.F116.1,by="SNP.T")
FarmCPUCV0.05.F116.flo.S <- separate(FarmCPUCV0.05.F116.flo, SNP.T, c("SNP","Trait.name"),sep="/")

str(FarmCPUCV0.05.F116.flo.S)


###for flowering 96
###for flowering 96
###for flowering 96
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYf.96.4814.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDf.96.4814.csv")
flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
flo <- flo.1[flo.1$Ind.SNP=="F96+4814" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.F96")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.F96/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.F96 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.F96$filename)

FarmCPUCV0.05.F96$filename <- gsub('FarmCPUCV0.05.F96/Trait.', '', FarmCPUCV0.05.F96$filename)
FarmCPUCV0.05.F96$filename <- gsub('.csv', '', FarmCPUCV0.05.F96$filename)

levels(as.factor(FarmCPUCV0.05.F96$filename))

FarmCPUCV0.05.F96.1 <- unite(FarmCPUCV0.05.F96,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.F96.flo <- merge(flo.2,FarmCPUCV0.05.F96.1,by="SNP.T")
FarmCPUCV0.05.F96.flo.S <- separate(FarmCPUCV0.05.F96.flo, SNP.T, c("SNP","Trait.name"),sep="/")


###for Culm 106
###for Culm 106
###for Culm 106
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.106.4202.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.106.4202.csv")
flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C106+4202" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.C106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.C106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.C106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.C106$filename)

FarmCPUCV0.05.C106$filename <- gsub('FarmCPUCV0.05.C106/Trait.', '', FarmCPUCV0.05.C106$filename)
FarmCPUCV0.05.C106$filename <- gsub('.csv', '', FarmCPUCV0.05.C106$filename)

levels(as.factor(FarmCPUCV0.05.C106$filename))

FarmCPUCV0.05.C106.1 <- unite(FarmCPUCV0.05.C106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.C106.flo <- merge(flo.2,FarmCPUCV0.05.C106.1,by="SNP.T")
FarmCPUCV0.05.C106.flo.S <- separate(FarmCPUCV0.05.C106.flo, SNP.T, c("SNP","Trait.name"),sep="/")


###for Culm 116
###for Culm 116
###for Culm 116
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.116.3293.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.116.3293.csv")
flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C116+3293" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.C116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.C116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.C116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.C116$filename)

FarmCPUCV0.05.C116$filename <- gsub('FarmCPUCV0.05.C116/Trait.', '', FarmCPUCV0.05.C116$filename)
FarmCPUCV0.05.C116$filename <- gsub('.csv', '', FarmCPUCV0.05.C116$filename)

levels(as.factor(FarmCPUCV0.05.C116$filename))

FarmCPUCV0.05.C116.1 <- unite(FarmCPUCV0.05.C116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.C116.flo <- merge(flo.2,FarmCPUCV0.05.C116.1,by="SNP.T")
FarmCPUCV0.05.C116.flo.S <- separate(FarmCPUCV0.05.C116.flo, SNP.T, c("SNP","Trait.name"),sep="/")

###for Culm 124
###for Culm 124
###for Culm 124
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.124.2560.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.124.2560.csv")
flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C124+2560" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.C124")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.C124/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.C124 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.C124$filename)

FarmCPUCV0.05.C124$filename <- gsub('FarmCPUCV0.05.C124/Trait.', '', FarmCPUCV0.05.C124$filename)
FarmCPUCV0.05.C124$filename <- gsub('.csv', '', FarmCPUCV0.05.C124$filename)

levels(as.factor(FarmCPUCV0.05.C124$filename))

FarmCPUCV0.05.C124.1 <- unite(FarmCPUCV0.05.C124,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.C124.flo <- merge(flo.2,FarmCPUCV0.05.C124.1,by="SNP.T")
FarmCPUCV0.05.C124.flo.S <- separate(FarmCPUCV0.05.C124.flo, SNP.T, c("SNP","Trait.name"),sep="/")


##for OWA 102
###for OWA 102
###for OWA 102
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.102.4322.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.102.4322.csv")


flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")

flo <- flo.1[flo.1$Ind.SNP=="O102+4322" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Entry"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.O102")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.O102/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

FarmCPUCV0.05.O102 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.O102$filename)

FarmCPUCV0.05.O102$filename <- gsub('FarmCPUCV0.05.O102/Trait.', '', FarmCPUCV0.05.O102$filename)
FarmCPUCV0.05.O102$filename <- gsub('.csv', '', FarmCPUCV0.05.O102$filename)

levels(as.factor(FarmCPUCV0.05.O102$filename))

FarmCPUCV0.05.O102.1 <- unite(FarmCPUCV0.05.O102,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.O102.flo <- merge(flo.2,FarmCPUCV0.05.O102.1,by="SNP.T")
FarmCPUCV0.05.O102.flo.S <- separate(FarmCPUCV0.05.O102.flo, SNP.T, c("SNP","Trait.name"),sep="/")



##for OWA 112
###for OWA 112
###for OWA 112
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.112.3450.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.112.3450.csv")
flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")

flo <- flo.1[flo.1$Ind.SNP=="O112+3450" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Entry"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.O112")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.O112/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

FarmCPUCV0.05.O112 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.O112$filename)

FarmCPUCV0.05.O112$filename <- gsub('FarmCPUCV0.05.O112/Trait.', '', FarmCPUCV0.05.O112$filename)
FarmCPUCV0.05.O112$filename <- gsub('.csv', '', FarmCPUCV0.05.O112$filename)

levels(as.factor(FarmCPUCV0.05.O112$filename))

FarmCPUCV0.05.O112.1 <- unite(FarmCPUCV0.05.O112,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.O112.flo <- merge(flo.2,FarmCPUCV0.05.O112.1,by="SNP.T")
FarmCPUCV0.05.O112.flo.S <- separate(FarmCPUCV0.05.O112.flo, SNP.T, c("SNP","Trait.name"),sep="/")


##for OWA 122
###for OWA 122
###for OWA 122
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.122.2646.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.122.2646.csv")
flo.1 <- read.csv("GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")

flo <- flo.1[flo.1$Ind.SNP=="O122+2646" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Entry"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.O122")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.O122/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}
#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.O122 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.O122$filename)

FarmCPUCV0.05.O122$filename <- gsub('FarmCPUCV0.05.O122/Trait.', '', FarmCPUCV0.05.O122$filename)
FarmCPUCV0.05.O122$filename <- gsub('.csv', '', FarmCPUCV0.05.O122$filename)

levels(as.factor(FarmCPUCV0.05.O122$filename))

FarmCPUCV0.05.O122.1 <- unite(FarmCPUCV0.05.O122,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.O122.flo <- merge(flo.2,FarmCPUCV0.05.O122.1,by="SNP.T")
FarmCPUCV0.05.O122.flo.S <- separate(FarmCPUCV0.05.O122.flo, SNP.T, c("SNP","Trait.name"),sep="/")



###combine all of the result
FarmCPUCV0.05.PVE.a <- do.call("rbind",list(FarmCPUCVim.flo.S,FarmCPUCV0.05.C106.flo.S,FarmCPUCV0.05.C116.flo.S,FarmCPUCV0.05.C124.flo.S,
                                        FarmCPUCV0.05.F96.flo.S,FarmCPUCV0.05.O102.flo.S,
                                        FarmCPUCV0.05.O112.flo.S,FarmCPUCV0.05.O122.flo.S))

write.csv(FarmCPUCV0.05.PVE.a,file="GAPIT0.05Result/FarmCPUCV0.05.PVE.a.1.csv")



