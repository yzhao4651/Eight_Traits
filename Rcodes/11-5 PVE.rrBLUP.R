### this one for MLMSUPER
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
flo.1 <- read.csv("GAPIT0.05Result/rrBLUPalltrait.PCA.csv")
str(flo.1)
levels(flo.1$Ind.SNP)
#colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
levels(flo.1$Trait.name)
flo <- flo.1[flo.1$Ind.SNP=="124im+36088",]

flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.15.116im")
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.25.116im")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.36088imup")
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
fileNames <- Sys.glob(paste("rrBLUPPCA.36088imup/*.", extension, sep = ""))

#myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")

mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCAim <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCAim$filename)

rrBLUPPCAim$filename <- gsub('rrBLUPPCA.36088imup/Trait.', '', rrBLUPPCAim$filename)
rrBLUPPCAim$filename <- gsub('.csv', '', rrBLUPPCAim$filename)

levels(as.factor(rrBLUPPCAim$filename))

rrBLUPPCAim <- unite(rrBLUPPCAim,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCAim.flo <- merge(flo,rrBLUPPCAim,by="SNP.T")
rrBLUPPCAim.flo.S <- separate(rrBLUPPCAim.flo, SNP.T, c("SNP","Trait.name"),sep="/")
str(rrBLUPPCAim.flo.S)


###for flowering 106
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
flo.1 <- read.csv("GAPIT0.05Result/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="F106+4185",]
flo <- droplevels(flo)
levels(flo.$Trait.name)
flo <- flo[!(flo$Trait.name=="Surv"),]
levels(flo$Trait.name)

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.F106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.F106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.F106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.F106$filename)

rrBLUPPCA.F106$filename <- gsub('rrBLUPPCA.F106/Trait.', '', rrBLUPPCA.F106$filename)
rrBLUPPCA.F106$filename <- gsub('.csv', '', rrBLUPPCA.F106$filename)

levels(as.factor(rrBLUPPCA.F106$filename))

rrBLUPPCA.F106 <- unite(rrBLUPPCA.F106,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.F106.flo <- merge(flo,rrBLUPPCA.F106,by="SNP.T")
rrBLUPPCA.F106.flo.S <- separate(rrBLUPPCA.F106.flo, SNP.T, c("SNP","Trait.name"),sep="/")


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
flo.1 <- read.csv("GAPIT0.05Result/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="C106+4202",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.C106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.C106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.C106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.C106$filename)

rrBLUPPCA.C106$filename <- gsub('rrBLUPPCA.C106/Trait.', '', rrBLUPPCA.C106$filename)
rrBLUPPCA.C106$filename <- gsub('.csv', '', rrBLUPPCA.C106$filename)

levels(as.factor(rrBLUPPCA.C106$filename))

rrBLUPPCA.C106 <- unite(rrBLUPPCA.C106,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.C106.flo <- merge(flo,rrBLUPPCA.C106,by="SNP.T")
rrBLUPPCA.C106.flo.S <- separate(rrBLUPPCA.C106.flo, SNP.T, c("SNP","Trait.name"),sep="/")


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
flo.1 <- read.csv("GAPIT0.05Result/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="C116+3293",]
str(flo)
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.C116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}



setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.C116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  str(sample)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.C116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.C116$filename)

rrBLUPPCA.C116$filename <- gsub('rrBLUPPCA.C116/Trait.', '', rrBLUPPCA.C116$filename)
rrBLUPPCA.C116$filename <- gsub('.csv', '', rrBLUPPCA.C116$filename)

levels(as.factor(rrBLUPPCA.C116$filename))

rrBLUPPCA.C116 <- unite(rrBLUPPCA.C116,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.C116.flo <- merge(flo,rrBLUPPCA.C116,by="SNP.T")
rrBLUPPCA.C116.flo.S <- separate(rrBLUPPCA.C116.flo, SNP.T, c("SNP","Trait.name"),sep="/")

###for Culm 124
###for Culm 124
###for Culm 124
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.124.2560.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.124.2560.csv")
flo.1 <- read.csv("GAPIT0.05Result/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="C124+2560",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.C124")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.C124/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.C124 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.C124$filename)

rrBLUPPCA.C124$filename <- gsub('rrBLUPPCA.C124/Trait.', '', rrBLUPPCA.C124$filename)
rrBLUPPCA.C124$filename <- gsub('.csv', '', rrBLUPPCA.C124$filename)

levels(as.factor(rrBLUPPCA.C124$filename))

rrBLUPPCA.C124 <- unite(rrBLUPPCA.C124,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.C124.flo <- merge(flo,rrBLUPPCA.C124,by="SNP.T")
rrBLUPPCA.C124.flo.S <- separate(rrBLUPPCA.C124.flo, SNP.T, c("SNP","Trait.name"),sep="/")

####

###for O102
###for O102
###for O102 O102+4322
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.102.4322.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.102.4322.csv")
flo.1 <- read.csv("GAPIT0.05Result/rrBLUPalltrait.PCA.csv")
levels(flo.1$Trait.name)
levels(flo.1$Ind.SNP)
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="O102+4322",]
flo <- droplevels(flo)
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.O102")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.O102/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.O102 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.O102$filename)

rrBLUPPCA.O102$filename <- gsub('rrBLUPPCA.O102/Trait.', '', rrBLUPPCA.O102$filename)
rrBLUPPCA.O102$filename <- gsub('.csv', '', rrBLUPPCA.O102$filename)

levels(as.factor(rrBLUPPCA.O102$filename))

rrBLUPPCA.O102 <- unite(rrBLUPPCA.O102,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.O102.flo <- merge(flo,rrBLUPPCA.O102,by="SNP.T")
rrBLUPPCA.O102.flo.S <- separate(rrBLUPPCA.O102.flo, SNP.T, c("SNP","Trait.name"),sep="/")

###for O112
###for O112
###for O112 O112+3450
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.112.3450.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.112.3450.csv")
flo.1 <- read.csv("GAPIT0.05Result/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="O112+3450",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.O112")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.O112/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.O112 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.O112$filename)

rrBLUPPCA.O112$filename <- gsub('rrBLUPPCA.O112/Trait.', '', rrBLUPPCA.O112$filename)
rrBLUPPCA.O112$filename <- gsub('.csv', '', rrBLUPPCA.O112$filename)

levels(as.factor(rrBLUPPCA.O112$filename))

rrBLUPPCA.O112 <- unite(rrBLUPPCA.O112,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.O112.flo <- merge(flo,rrBLUPPCA.O112,by="SNP.T")
rrBLUPPCA.O112.flo.S <- separate(rrBLUPPCA.O112.flo, SNP.T, c("SNP","Trait.name"),sep="/")


####O122
###for O122
###for O122
###for O122 O122+2646
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.122.2646.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.122.2646.csv")
flo.1 <- read.csv("GAPIT0.05Result/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="O122+2646",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.O122")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.O122/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.O122 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.O122$filename)

rrBLUPPCA.O122$filename <- gsub('rrBLUPPCA.O122/Trait.', '', rrBLUPPCA.O122$filename)
rrBLUPPCA.O122$filename <- gsub('.csv', '', rrBLUPPCA.O122$filename)

levels(as.factor(rrBLUPPCA.O122$filename))

rrBLUPPCA.O122 <- unite(rrBLUPPCA.O122,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.O122.flo <- merge(flo,rrBLUPPCA.O122,by="SNP.T")
rrBLUPPCA.O122.flo.S <- separate(rrBLUPPCA.O122.flo, SNP.T, c("SNP","Trait.name"),sep="/")
str(rrBLUPPCA.O122.flo.S)

rrBLUPPCA.O122.flo.S$PVE <- as.character(rrBLUPPCA.O122.flo.S$PVE)
rrBLUPPCA.O112.flo.S$PVE <- as.character(rrBLUPPCA.O112.flo.S$PVE)

###combine all of them together 

rrBLUPPCA.PVE.a <- do.call("rbind",list(rrBLUPPCA.O102.flo.S, rrBLUPPCA.O112.flo.S,rrBLUPPCA.O122.flo.S,rrBLUPPCA.C106.flo.S,
                                        rrBLUPPCA.C116.flo.S,rrBLUPPCA.C124.flo.S,
                                        rrBLUPPCA.F106.flo.S,rrBLUPPCAim.flo.S))

write.csv(rrBLUPPCA.PVE.a,file="GAPIT0.05Result/rrBLUPPCA.PVE.a.csv")



