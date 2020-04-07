
###import all the dataset with MAF 0.05
###ECMML method with CV
EPS <- read.csv("GAPIT0.05Result/ECMML0.05.PS.csv",row.names=1)
str(EPS)
#chang maf to MAF
colnames(EPS)[colnames(EPS)=="maf"] <- "MAF"
str(EPS)
#EPSNO <- read.csv("GAPIT0.05Result/ECMML0.05.Adj.P.PSNO.csv",row.names=1)
#colnames(EPSNO)[colnames(EPSNO)=="maf"] <- "MAF"
#str(EPSNO)
###method with GLM, MLM, MLMM with CV
GPS <- read.csv("GAPIT0.05Result/GLMMLMMLMM0.05.PS.csv",row.names=1)
colnames(GPS)[colnames(GPS)=="maf"] <- "MAF"
str(GPS)
#####method SUPER
CMLM<- read.csv("GAPIT0.05Result/CMLMSUPER0.05.PS.PVE.2.csv",row.names=1)
colnames(CMLM)[colnames(CMLM)=="maf"] <- "MAF"
str(CMLM)
#####method SUPER
MLM <- read.csv("GAPIT0.05Result/MLMSUPER0.05.PS.PVE.a.2.csv",row.names=1)
colnames(MLM)[colnames(MLM)=="maf"] <- "MAF"
str(MLM)

####rMLMM with flowering with CV
rMFPS <- read.csv("GAPIT0.05Result/rMLMMflo.PS.csv",row.names=1)
str(rMFPS)
####rMLMM with flowering with NO CV
rMFPSNO <- read.csv("GAPIT0.05Result/rMLMMflo.PSNO.csv",row.names=1)
str(rMFPSNO)
###mrMLM with Culm trait with CV
rMCPS <- read.csv("GAPIT0.05Result/rMLMMCulm.PS.csv",row.names=1)
str(rMFPS)
###mrMLM with Culm trait with NO CV
rMCPSNO <- read.csv("GAPIT0.05Result/rMLMMCulm.PSNO.csv",row.names=1)
str(rMFPSNO)
rrBLUP <- read.csv("GAPIT0.05Result/rrBLUPPCA.PVE.a.csv",row.names=1)
str(rrBLUP)
###FarmCPU without CV
FarmPSNO <- read.csv("GAPIT0.05Result/FarmCPU0.05.PVE.a.1.csv",row.names=1)
colnames(FarmPSNO)[colnames(FarmPSNO)=="maf"] <- "MAF"
str(FarmPSNO)
####FarmCPU without CV
FarmPS <- read.csv("GAPIT0.05Result/FarmCPUCV0.05.PVE.a.1.csv",row.names=1)
colnames(FarmPS)[colnames(FarmPS)=="maf"] <- "MAF"
str(FarmPS)


###combineing all of them together 
alltrait.method <- plyr::ldply(list(EPS,GPS,rMFPS,rMFPSNO,rMCPS,rMCPSNO,rrBLUP,FarmPS,FarmPSNO,CMLM,MLM))
str(alltrait.method)

###
###change the P.value to with zero decime 
alltrait.method$P.value.1 <- as.numeric(formatC(alltrait.method$P.value, format = "e", digits = 0))
#alltrait.method$MAF.1<- as.numeric(formatC(alltrait.method$MAF, format = "e", digits = 3))
#alltrait.method$PVE.1 <- as.numeric(formatC(alltrait.method$PVE, format = "e", digits = 4))
str(alltrait.method)
##write out the total traits 
write.csv(alltrait.method, file="data/alltraitsGWASresultupdated1.csv")

###select all the flowering trait

###import the all traits 
all <- read.csv("data/alltraitsGWASresultupdated1.csv")
levels(all$Trait.name)
all <- alltrait.method
floall <- all[all$Trait.name =="HD_1" | all$Trait.name =="HD_50" | 
                               all$Trait.name =="FD_1" | all$Trait.name =="FD_50" |
                               all$Trait.name =="FM_1" | all$Trait.name =="FM_50" |
                               all$Trait.name =="HM_1" | all$Trait.name =="HM_50" |
                               all$Trait.name =="FW_1" | all$Trait.name =="FW_50" |
                               all$Trait.name =="HW_1" | all$Trait.name =="HW_50" |
                               all$Trait.name =="GHW_1" | all$Trait.name =="GHW_50" |
                               all$Trait.name =="GFW_1" | all$Trait.name =="GFW_50" |
                               all$Trait.name =="fprind" | all$Trait.name =="fprinGW" |
                               all$Trait.name =="fprinM" | all$Trait.name =="fprinW",]

floall <- droplevels(floall)
###write out the flowering trait 
levels(floall$Ind.SNP)
floall$Ind.SNP <- as.factor(floall$Ind.SNP)
levels(floall$Ind.SNP)[levels(floall$Ind.SNP)=="124+36088im"] <- "F116+36088im"
levels(floall$Ind.SNP)[levels(floall$Ind.SNP)=="124im+36088"] <- "F116+36088im"
levels(floall$Ind.SNP)[levels(floall$Ind.SNP)=="F116+36088im"] <- "F116+36088im"
#floall$Ind.SNP <- as.factor(floall$Ind.SNP)
levels(floall$Ind.SNP)
str(floall)
floall.1 <- floall[,c(1:5,12:17,24)]
write.csv(floall.1,file="Flotraits/floalltraits2.csv")
###find the largest PVE and MAF

### select the largest MAF and PVE
floall <- read.csv("Flotraits/floalltraits2.csv",row.names = 1)
str(floall)
floall <- floall[order(floall$SNP, -(as.numeric(floall$MAF))),] #sort by id and reverse of abs(value)
floall.MAF.2 <- floall[ !duplicated(floall$SNP), ] 
floall <- read.csv("Flotraits/floalltraits1.csv",row.names = 1)
str(floall) 
floall <- floall[order(floall$SNP, -(as.numeric(floall$PVE))), ] #sort by id and reverse of abs(value)
floall.PVE.2 <- floall[ !duplicated(floall$SNP), ] 
##merge this two 
flo.MAF.PVE.Max <- merge(x=floall.MAF.2,y=floall.PVE.2, by="SNP" )
str(flo.MAF.PVE.Max)

flo.MAF.PVE.Max.1 <- flo.MAF.PVE.Max[,c(1,5,34)]

flo.MAF.PVE.Max.1$ MAF.x <- as.numeric(formatC(flo.MAF.PVE.Max.1$ MAF.x, format = "e", digits = 3))
flo.MAF.PVE.Max.1$ PVE.y <- as.numeric(formatC(flo.MAF.PVE.Max.1$ PVE.y, format = "e", digits = 3))

write.csv(flo.MAF.PVE.Max.1,file="Flotraits/flo.MAF.PVE.Max.2.csv")



###select the flo
##import all of the data 
all <- read.csv("data/alltraitsGWASresultupdated1.csv")
all <- alltrait.method 
levels(all$Trait.name)
levels(all$Trait.name)[levels(all$Trait.name)=="OWA.18"] <- "OWA18"
levels(all$Trait.name)[levels(all$Trait.name)=="OWA.19"] <- "OWA19"
levels(all$Ind.SNP)[levels(all$Ind.SNP)=="124im+36088"] <- "122+35256im"
levels(all$Ind.SNP)[levels(all$Ind.SNP)=="124+36088im"] <- "122+35256im"
OWAall <- all[all$Trait.name =="OWA18"|all$Trait.name =="OWA19",] 
str(OWAall)
OWAall <- droplevels(OWAall)
##write out the owa data set 

write.csv(OWAall,file="OWA2/OWA2.csv") 

### select the largest MAF and PVE

OWAall <- OWAall[order(OWAall$SNP, -(as.numeric(OWAall$MAF))),] #sort by id and reverse of abs(value)
OWAall.MAF.2 <- OWAall[ !duplicated(OWAall$SNP), ] 
OWAall <- OWAall[order(OWAall$SNP, -(as.numeric(OWAall$PVE))), ] #sort by id and reverse of abs(value)
OWAall.PVE.2 <- OWAall[ !duplicated(OWAall$SNP), ] 
##merge this two 
OWA.MAF.PVE.Max <- merge(x=OWAall.MAF.2,y=OWAall.PVE.2, by="SNP" )
str(OWA.MAF.PVE.Max)
OWA.MAF.PVE.Max.1 <- OWA.MAF.PVE.Max[,c(1,5,40)]
str(OWA.MAF.PVE.Max.1)

#OWA.MAF.PVE.Max.1$ MAF.x <- as.numeric(formatC(OWA.MAF.PVE.Max.1$ MAF.x, format = "e", digits = 3))
#OWA.MAF.PVE.Max.1$ PVE.y <- as.numeric(formatC(OWA.MAF.PVE.Max.1$ PVE.y, format = "e", digits = 3))

write.csv(OWA.MAF.PVE.Max.1,file="OWA2/OWA2.MAF.PVE.Max.1.csv")

OWAall <- all[all$Trait.name =="OWA18"|all$Trait.name =="OWA19",] 
str(OWAall)
OWAall <- droplevels(OWAall)
OWA.MAF.PVE.Max <- list()
for (i in levels(OWAall$Trait.name)){
  OWAall.1 <- OWAall[OWAall$Trait.name==i,]
  OWAall.2 <- OWAall.1[order(OWAall.1$SNP, -(as.numeric(OWAall.1$MAF))),] #sort by id and reverse of abs(value)
  OWAall.MAF.2 <- OWAall.2[ !duplicated(OWAall.2$SNP), ] 
  OWAall.3 <- OWAall.1[order(OWAall.1$SNP, -(as.numeric(OWAall.1$PVE))), ] #sort by id and reverse of abs(value)
  OWAall.PVE.2 <- OWAall.3[ !duplicated(OWAall.3$SNP), ] 
  OWA.MAF.PVE.Max [[i]]<- merge(x=OWAall.MAF.2,y=OWAall.PVE.2, by="SNP" )
}
OWA.MAF.PVE.Max <- do.call(rbind,OWA.MAF.PVE.Max)
str(OWA.MAF.PVE.Max)

####write the data out set 
write.csv(OWA.MAF.PVE.Max,file="OWA2/OWA2.MAF.PVE.Max.2.csv")

###import the all traits 
all <- read.csv("GAPIT0.05Result/Alltrait.1-org.1.csv")
levels(all$Trait.name)
###get SRD traits 
SRDall <- all[all$Trait.name =="SRD",] 
##write out the owa data set 
write.csv(SRDall,file="OWA/SRD.csv")
##gwet the large number of MAF 

SRDall <- SRDall[order(SRDall$SNP, -(as.numeric(SRDall$MAF.1))),] #sort by id and reverse of abs(value)
SRDall.MAF.2 <- SRDall[ !duplicated(SRDall$SNP), ] 
SRDall <- SRDall[order(SRDall$SNP, -(as.numeric(SRDall$PVE.1))), ] #sort by id and reverse of abs(value)
SRDall.PVE.2 <- SRDall[ !duplicated(SRDall$SNP), ] 
##merge this two 
write.csv(SRDall.PVE.2,file="OWA/SRD.MAF.csv")
###get all the SNP for trait ADD
all <- read.csv("GAPIT0.05Result/Alltrait.1-org.1.csv")
levels(all$Trait.name)
###get ADD traits 
ADDall <- all[all$Trait.name =="ADD",]
##write out the owa data set 
write.csv(ADDall,file="OWA/ADD.csv")

##get the large number of MAF 
ADDall <- ADDall[order(ADDall$SNP, -(as.numeric(ADDall$MAF.1))),] #sort by id and reverse of abs(value)
ADDall.MAF.2 <- ADDall[ !duplicated(ADDall$SNP), ] 
ADDall <- ADDall[order(ADDall$SNP, -(as.numeric(ADDall$PVE.1))), ] #sort by id and reverse of abs(value)
ADDall.PVE.2 <- ADDall[ !duplicated(ADDall$SNP), ] 
write.csv(ADDall.PVE.2,file="OWA/ADD.MAF.csv")



  
  
  
  
  
  
  
  
  
  