

###import the SRD and OWA 
myYSRD <- read.csv("data/myYimputedSNP19SRD.csv")
myYOWA <- read.csv("data/myYimputedSNP19OWA2.csv")
myY1 <- read.csv("data/myYimputedSNP19.csv")
names(myY1)
myY1 <- myY1[,-c(2,18)]
myY <- Reduce(function(x, y) merge(x, y, all.x=TRUE), list(myY1,myYOWA,myYSRD))
names(myY)
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
flo.1 <- read.csv("GAPIT0.05Result/SUPER0.05.PS.csv")
str(flo.1)

flo <- flo.1[flo.1$Ind.SNP=="124+36088im" & flo.1$Method=="CMLM+SUPER",]
str(flo)
names(flo)
levels(flo$Trait.name)
flo <- droplevels(flo)
levels(flo$Ind.SNP)
levels(as.factor(flo$Method))
levels(as.factor(flo$Trait.name))
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="OWA"),]
flo <- flo[!(flo$Trait.name=="fprin"),]


#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.15.116im")
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.8.3.116im")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.36088imUP2")
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
fileNames <- Sys.glob(paste("CMLMSUPER0.05.36088imUP2/*.", extension, sep = ""))

#myY <- read.csv("data/myYimputedSNP19.csv")

myGD <- read.csv("data/myGDimputedSNP19.csv")

mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMCUPERim <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMCUPERim$filename)

CMLMCUPERim$filename <- gsub('CMLMSUPER0.05.36088imUP2/Trait.', '', CMLMCUPERim$filename)
CMLMCUPERim$filename <- gsub('.csv', '', CMLMCUPERim$filename)

levels(as.factor(CMLMCUPERim$filename))
str(CMLMCUPERim)
CMLMCUPERim.1 <- unite(CMLMCUPERim,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

CMLMCUPERim.flo <- merge(flo.2,CMLMCUPERim.1,by="SNP.T")

CMLMCUPERim.flo.S <- separate(CMLMCUPERim.flo, SNP.T, c("SNP","Trait.name"),sep="/")
str(CMLMCUPERim.flo.S)

write.csv(CMLMCUPERim.flo.S,file="GAPIT0.05Result/CMLMCUPERim.flo.S2.csv")



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

flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
str(flo.1)

flo <- flo.1[flo.1$Ind.SNP=="F106+4185" & flo.1$Method=="CMLM+SUPER",]

flo <- flo[!(flo$Trait.name=="Surv"),]

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.F106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("CMLMSUPER0.05.F106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMSUPER0.05.F106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMSUPER0.05.F106$filename)

CMLMSUPER0.05.F106$filename <- gsub('CMLMSUPER0.05.F106/Trait.', '', CMLMSUPER0.05.F106$filename)
CMLMSUPER0.05.F106$filename <- gsub('.csv', '', CMLMSUPER0.05.F106$filename)

levels(as.factor(CMLMSUPER0.05.F106$filename))

CMLMSUPER0.05.F106.1 <- unite(CMLMSUPER0.05.F106,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

CMLMSUPER0.05.F106.flo <- merge(flo.2,CMLMSUPER0.05.F106.1,by="SNP.T")
CMLMSUPER0.05.F106.flo.S <- separate(CMLMSUPER0.05.F106.flo, SNP.T, c("SNP","Trait.name"),sep="/")
CMLMSUPER0.05.F106.flo.S <- CMLMSUPER0.05.F106.flo.S[-18]
str(CMLMSUPER0.05.F106.flo.S)
write.csv(CMLMCUPERim.flo.S,file="GAPIT0.05Result/CMLMSUPER0.05.F106.flo.S.csv")

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

flo.1 <- read.csv("GAPIT0.05Result/SUPER0.05.PS.csv")
str(flo.1)

flo <- flo.1[flo.1$Ind.SNP=="F116+3098" & flo.1$Method=="CMLM+SUPER",]

flo <- flo[!(flo$Trait.name=="Surv"),]

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.F116UP")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("CMLMSUPER0.05.F116UP/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMSUPER0.05.F116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMSUPER0.05.F116$filename)

CMLMSUPER0.05.F116$filename <- gsub('CMLMSUPER0.05.F116UP/Trait.', '', CMLMSUPER0.05.F116$filename)
CMLMSUPER0.05.F116$filename <- gsub('.csv', '', CMLMSUPER0.05.F116$filename)

levels(as.factor(CMLMSUPER0.05.F116$filename))

CMLMSUPER0.05.F116.1 <- unite(CMLMSUPER0.05.F116,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

CMLMSUPER0.05.F116.flo <- merge(flo.2,CMLMSUPER0.05.F116.1,by="SNP.T")
CMLMSUPER0.05.F116.flo.S <- separate(CMLMSUPER0.05.F116.flo, SNP.T, c("SNP","Trait.name"),sep="/")
write.csv(CMLMCUPERim.flo.S,file="GAPIT0.05Result/CMLMSUPER0.05.F116.flo.S.csv")


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

flo.1 <- read.csv("GAPIT0.05Result/SUPER0.05.PS.csv")

str(flo.1)

flo <- flo.1[flo.1$Ind.SNP=="F96+4814" & flo.1$Method=="CMLM+SUPER",]
str(flo)

flo <- flo[!(flo$Trait.name=="Surv"),]

str(flo)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.F96")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("CMLMSUPER0.05.F96/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMSUPER0.05.F96 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMSUPER0.05.F96$filename)

CMLMSUPER0.05.F96$filename <- gsub('CMLMSUPER0.05.F96/Trait.', '', CMLMSUPER0.05.F96$filename)
CMLMSUPER0.05.F96$filename <- gsub('.csv', '', CMLMSUPER0.05.F96$filename)

levels(as.factor(CMLMSUPER0.05.F96$filename))

CMLMSUPER0.05.F96.1 <- unite(CMLMSUPER0.05.F96,"SNP.T", SNP, filename, sep="/")
flo.2 <- unite(flo,"SNP.T", SNP, Trait.name, sep="/")

CMLMSUPER0.05.F96.flo <- merge(flo.2,CMLMSUPER0.05.F96.1, by="SNP.T")
CMLMSUPER0.05.F96.flo.S <- separate(CMLMSUPER0.05.F96.flo, SNP.T, c("SNP","Trait.name"),sep="/")
str(CMLMSUPER0.05.F96.flo.S)

write.csv(CMLMCUPERim.flo.S,file="GAPIT0.05Result/CMLMSUPER0.05.F96.flo.S.csv")


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
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C106+4202" & flo.1$Method=="CMLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]


setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.C106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("CMLMSUPER0.05.C106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMSUPER0.05.C106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMSUPER0.05.C106$filename)

CMLMSUPER0.05.C106$filename <- gsub('CMLMSUPER0.05.C106/Trait.', '', CMLMSUPER0.05.C106$filename)
CMLMSUPER0.05.C106$filename <- gsub('.csv', '', CMLMSUPER0.05.C106$filename)

levels(as.factor(CMLMSUPER0.05.C106$filename))

CMLMSUPER0.05.C106.1 <- unite(CMLMSUPER0.05.C106,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

CMLMSUPER0.05.C106.flo <- merge(flo.2,CMLMSUPER0.05.C106.1,by="SNP.T")
CMLMSUPER0.05.C106.flo.S <- separate(CMLMSUPER0.05.C106.flo, SNP.T, c("SNP","Trait.name"),sep="/")
CMLMSUPER0.05.C106.flo.S <- CMLMSUPER0.05.C106.flo.S[-18]
write.csv(CMLMSUPER0.05.C106.flo.S,file="GAPIT0.05Result/CMLMSUPER0.05.C106.flo.S.csv")

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
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C116+3293" & flo.1$Method=="CMLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]


setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.C116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("CMLMSUPER0.05.C116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMSUPER0.05.C116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMSUPER0.05.C116$filename)

CMLMSUPER0.05.C116$filename <- gsub('CMLMSUPER0.05.C116/Trait.', '', CMLMSUPER0.05.C116$filename)
CMLMSUPER0.05.C116$filename <- gsub('.csv', '', CMLMSUPER0.05.C116$filename)

levels(as.factor(CMLMSUPER0.05.C116$filename))

CMLMSUPER0.05.C116.1 <- unite(CMLMSUPER0.05.C116,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

CMLMSUPER0.05.C116.flo.2 <- merge(flo.2,CMLMSUPER0.05.C116.1,by="SNP.T")
CMLMSUPER0.05.C116.flo.S <- separate(CMLMSUPER0.05.C116.flo.2, SNP.T, c("SNP","Trait.name"),sep="/")
CMLMSUPER0.05.C116.flo.S <- CMLMSUPER0.05.C116.flo.S[-18]
str(CMLMSUPER0.05.C116.flo.S)
write.csv(CMLMSUPER0.05.C116.flo.S,file="GAPIT0.05Result/CMLMSUPER0.05.C116.flo.S.csv")

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

#flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo.1 <- read.csv("GAPIT0.05Result/SUPER0.05.PS.csv")
str(flo.1)
flo <- flo.1[flo.1$Ind.SNP=="C124+2560" & flo.1$Method=="CMLM+SUPER",]
#flo <- flo.1[flo.1$Ind.SNP=="124.2560" & flo.1$Method=="CMLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]


setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.C124")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("CMLMSUPER0.05.C124/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMSUPER0.05.C124 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMSUPER0.05.C124$filename)

CMLMSUPER0.05.C124$filename <- gsub('CMLMSUPER0.05.C124/Trait.', '', CMLMSUPER0.05.C124$filename)
CMLMSUPER0.05.C124$filename <- gsub('.csv', '', CMLMSUPER0.05.C124$filename)

levels(as.factor(CMLMSUPER0.05.C124$filename))

CMLMSUPER0.05.C124.1 <- unite(CMLMSUPER0.05.C124,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

CMLMSUPER0.05.C124.flo <- merge(flo.2,CMLMSUPER0.05.C124.1,by="SNP.T")
CMLMSUPER0.05.C124.flo.S <- separate(CMLMSUPER0.05.C124.flo, SNP.T, c("SNP","Trait.name"),sep="/")
str(CMLMSUPER0.05.C116.flo.S)
write.csv(CMLMSUPER0.05.C116.flo.S,file="GAPIT0.05Result/CMLMSUPER0.05.C116.flo.S.csv")

###for OWA 102 this one no result 
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
flo.1 <- read.csv("GAPIT0.05Result/SUPER0.05.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="O102+4322" & flo.1$Method=="CMLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]


setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.O102")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("CMLMSUPER0.05.O102/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMSUPER0.05.O102 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMSUPER0.05.O102$filename)

CMLMSUPER0.05.O102$filename <- gsub('CMLMSUPER0.05.O102/Trait.', '', CMLMSUPER0.05.O102$filename)
CMLMSUPER0.05.O102$filename <- gsub('.csv', '', CMLMSUPER0.05.O102$filename)

levels(as.factor(CMLMSUPER0.05.O102$filename))

CMLMSUPER0.05.O102.1 <- unite(CMLMSUPER0.05.O102,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

CMLMSUPER0.05.O102.flo <- merge(flo.2,CMLMSUPER0.05.O102.1,by="SNP.T")
CMLMSUPER0.05.O102.flo.S <- separate(CMLMSUPER0.05.O102.flo, SNP.T, c("SNP","Trait.name"),sep="/")
write.csv(CMLMSUPER0.05.O102.flo.S,file="GAPIT0.05Result/CMLMSUPER0.05.O102.flo.S.csv")



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
flo.1 <- read.csv("GAPIT0.05Result/SUPER0.05.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="O112+3450" & flo.1$Method=="CMLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]


setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.O112")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("CMLMSUPER0.05.O112/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMSUPER0.05.O112 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMSUPER0.05.O112$filename)

CMLMSUPER0.05.O112$filename <- gsub('CMLMSUPER0.05.O112/Trait.', '', CMLMSUPER0.05.O112$filename)
CMLMSUPER0.05.O112$filename <- gsub('.csv', '', CMLMSUPER0.05.O112$filename)

levels(as.factor(CMLMSUPER0.05.O112$filename))

CMLMSUPER0.05.O112.1 <- unite(CMLMSUPER0.05.O112,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

CMLMSUPER0.05.O112.flo <- merge(flo.2,CMLMSUPER0.05.O112.1,by="SNP.T")
CMLMSUPER0.05.O112.flo.S <- separate(CMLMSUPER0.05.O112.flo, SNP.T, c("SNP","Trait.name"),sep="/")
write.csv(CMLMSUPER0.05.O112.flo.S,file="GAPIT0.05Result/CMLMSUPER0.05.O112.flo.S.csv")

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
flo.1 <- read.csv("GAPIT0.05Result/SUPER0.05.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="O122+2646" & flo.1$Method=="CMLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]


setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/CMLMSUPER0.05.O122")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("CMLMSUPER0.05.O122/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
CMLMSUPER0.05.O122 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(CMLMSUPER0.05.O122$filename)

CMLMSUPER0.05.O122$filename <- gsub('CMLMSUPER0.05.O122/Trait.', '', CMLMSUPER0.05.O122$filename)
CMLMSUPER0.05.O122$filename <- gsub('.csv', '', CMLMSUPER0.05.O122$filename)

levels(as.factor(CMLMSUPER0.05.O122$filename))

CMLMSUPER0.05.O122.1 <- unite(CMLMSUPER0.05.O122,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

CMLMSUPER0.05.O122.flo <- merge(flo.2,CMLMSUPER0.05.O122.1,by="SNP.T")
CMLMSUPER0.05.O122.flo.S <- separate(CMLMSUPER0.05.O122.flo, SNP.T, c("SNP","Trait.name"),sep="/")
str(CMLMSUPER0.05.O122.flo.S)
write.csv(CMLMSUPER0.05.O122.flo.S,file="GAPIT0.05Result/CMLMSUPER0.05.O122.flo.S.csv")

str(CMLMCUPERim.flo.S)
###combine all of file together 
CMLMSUPER0.05.PVE.a.2 <- do.call("rbind",list(CMLMSUPER0.05.C106.flo.S,CMLMSUPER0.05.C116.flo.S,CMLMSUPER0.05.C124.flo.S,
                                        CMLMSUPER0.05.O122.flo.S,CMLMSUPER0.05.O112.flo.S,CMLMCUPERim.flo.S,
                                        CMLMSUPER0.05.F106.flo.S,CMLMSUPER0.05.F116.flo.S,CMLMSUPER0.05.F96.flo.S))


write.csv(CMLMSUPER0.05.PVE.a.2,file="GAPIT0.05Result/CMLMSUPER0.05.PS.PVE.2.csv")

