###install the pcakages 
#install.packages("chromPlot")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("chromPlot")
library("chromPlot")

###import the GM as the gap dataset
load("~/Documents/MiacnathusSNPinformation /161025forGAPIT.RData")
##checking the GM 
str(myGM)
##changing the columan names of the this data set
colnames(myGM)[colnames(myGM)=="Chromosome"] <- "Chrom"
colnames(myGM)[colnames(myGM)=="Position"] <- "Start"
myGM$End <- myGM$Start
str(myGM)
myGM <- myGM[,c(2:4,1)]
str(myGM)


###draw the chromosomes with two groups 
##import the dataset from this study
##using all of the QTL and SNPs from this study and other study as bands.
###using the traits and who published as the ID to lable them on the chromosome

###for traits: Birce and BC and BBCC
###this one working 
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4 <- QTLID4[QTLID4$Group2=="BC"| QTLID4$Group2=="Bcirc"| QTLID4$Group2=="CCBC"|QTLID4$Group2=="CC/BC",]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")
#QTLID4 <- QTLID4[QTLID4$Group3=="Clark_et_al_(2019)" | QTLID4$Group3=="Current_Study",]
QTLID4 <- droplevels(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
str(QTLID4)
QTLID5 <- QTLID4[,c(1:3,6,11)]
str(QTLID5)
levels(as.factor(QTLID5$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
bands.bc <- bands.bc[bands.bc$Group2=="BC"| bands.bc$Group2=="Bcirc"| bands.bc$Group2=="CCBC"| bands.bc$Group2=="CC/BC",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc)
bands.bc <- droplevels(bands.bc)
bands.bc1 <- bands.bc[,c(1,4:5,10,12)]
bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc1)
levels(as.factor(bands.bc1$Chrom))
####using ID 
chromPlot(gaps=myGM, bands=bands.bc1, stat=QTLID5, statCol="Value",statName="Value")

###for traits: CCirc and CC and BBCC
###this one working 
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4.c <- QTLID4[QTLID4$Group2=="CC"| QTLID4$Group2=="CCirc"| QTLID4$Group2=="CCBC"|QTLID4$Group2=="CC/BC",]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")

#QTLID4.c <- QTLID4.c[QTLID4.c$Group3=="Clark_et_al_(2019)" | QTLID4.c$Group3=="Current_Study",]
QTLID4.c <- droplevels(QTLID4.c)
levels(QTLID4.c$Group2)
levels(QTLID4.c$Group3)
str(QTLID4.c)
QTLID5.c <- QTLID4.c[,c(1:3,6,11)]
str(QTLID5.c)
levels(as.factor(QTLID5.c$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
levels(bands.bc$Group2)
bands.bc1 <- bands.bc[bands.bc$Group2=="CC"| bands.bc$Group2=="CCirc"| bands.bc$Group2=="CCBC"| bands.bc$Group2=="CC/BC",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc)
bands.bc1 <- droplevels(bands.bc1)
bands.bc2 <- bands.bc1[,c(1,4:5,10,12)]
#bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc2)
levels(as.factor(bands.bc2$Chrom))
####using ID 
chromPlot(gaps=myGM, bands=bands.bc2, stat=QTLID5.c, statCol="Value",statName="Value")

####using the DBI
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1DBI.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4.BI <- QTLID4[QTLID4$Group2=="CmD_BI"| QTLID4$Group2=="DBI" | QTLID4$Group2=="stem_diameter"| QTLID4$Group2=="tiller_diameter"|
                      QTLID4$Group2=="CmVol"|QTLID4$Group2=="CmDW/V"|QTLID4$Group2=="CmV",]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")

#QTLID4.BI <- QTLID4.BI[QTLID4.BI$Group3=="Clark_et_al_(2019)" | QTLID4.BI$Group3=="Current_Study"| QTLID4.BI$Group3=="Clark_et_al_(2016)",]
QTLID4.BI <- droplevels(QTLID4.BI)
levels(QTLID4.BI$Group2)
levels(QTLID4.BI$Group3)
str(QTLID4.BI)
QTLID5.BI <- QTLID4.BI[,c(1:3,6,11)]
str(QTLID5.BI)
levels(as.factor(QTLID5.BI$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
levels(bands.bc$Group2)
bands.bc1 <- bands.bc[bands.bc$Group2=="CmD_BI"| bands.bc$Group2=="DBI"| bands.bc$Group2=="stem_diameter"|bands.bc$Group2=="tiller_diameter"|
                        bands.bc$Group2=="CmVol"| bands.bc$Group2=="CmDW/V"|bands.bc$Group2=="CmV",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc1)
bands.bc1 <- droplevels(bands.bc1)
levels(bands.bc1$Group3)
levels(bands.bc1$Group2)
bands.bc1 <- droplevels(bands.bc1)
bands.bc2 <- bands.bc1[,c(1,4:5,10,12)]
#bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc2)
levels(as.factor(bands.bc2$Chrom))
####using ID 

chromPlot(gaps=myGM, bands=bands.bc2, stat=QTLID5.BI, statCol="Value",statName="Value")

####using the DLI
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1DLI.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4.LI <- QTLID4[QTLID4$Group2=="CmD_LI"| QTLID4$Group2=="DTI" | QTLID4$Group2=="stem_diameter"| QTLID4$Group2=="tiller_diameter"| 
                      QTLID4$Group2=="CmV"|QTLID4$Group2=="CmDW/V" |QTLID4$Group2=="CmVol",,]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)

#QTLID4.LI <- QTLID4.LI[QTLID4.LI$Group3=="Clark_et_al_(2019)" | QTLID4.LI$Group3=="Current_Study" | QTLID4.BI$Group3=="Clark_et_al_(2016)",]
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)
str(QTLID4.LI)
QTLID5.LI <- QTLID4.LI[,c(1:3,6,11)]
str(QTLID5.LI)
levels(as.factor(QTLID5.LI$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
levels(bands.bc$Group2)
bands.bc1 <- bands.bc[bands.bc$Group2=="CmD_LI"| bands.bc$Group2=="DTI"| bands.bc$Group2=="stem_diameter"|bands.bc$Group2=="tiller_diameter"|
                        bands.bc$Group2=="CmV"| bands.bc$Group2=="CmDW/V"|bands.bc$Group2=="CmVol",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc)
bands.bc1 <- droplevels(bands.bc1)
bands.bc2 <- bands.bc1[,c(1,4:5,10,12)]
#bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc2)
levels(as.factor(bands.bc2$Chrom))
####using ID 

chromPlot(gaps=myGM, bands=bands.bc2, stat=QTLID5.LI, statCol="Value",statName="Value")

####using the Culm number
####using the CmN
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1CmN.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4.LI <- QTLID4[QTLID4$Group2=="CmN"| QTLID4$Group2=="IntL"|QTLID4$Group2=="CmNdN"|QTLID4$Group2=="CmNN"|QTLID4$Group2=="CmN", ]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)

#QTLID4.LI <- QTLID4.LI[QTLID4.LI$Group3=="Clark_et_al_(2019)" | QTLID4.LI$Group3=="Current_Study",]
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)
str(QTLID4.LI)
QTLID5.LI <- QTLID4.LI[,c(1:3,6,11)]
str(QTLID5.LI)
levels(as.factor(QTLID5.LI$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
levels(bands.bc$Group2)
bands.bc1 <- bands.bc[bands.bc$Group2=="CmN"| bands.bc$Group2=="IntL"|bands.bc$Group2=="CmNdN"|bands.bc$Group2=="CmNN"|bands.bc$Group2=="CmN",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc)
bands.bc1 <- droplevels(bands.bc1)
bands.bc2 <- bands.bc1[,c(1,4:5,10,12)]
#bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc2)
levels(as.factor(bands.bc2$Chrom))
####using ID 

chromPlot(gaps=myGM, bands=bands.bc2, stat=QTLID5.LI, statCol="Value",statName="Value")


####using the "Cml"
####
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1Cml.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4.LI <- QTLID4[QTLID4$Group2=="Cml"| QTLID4$Group2=="IntL"|QTLID4$Group2=="height"| QTLID4$Group2=="tallest_stem", ]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)

#QTLID4.LI <- QTLID4.LI[QTLID4.LI$Group3=="Clark_et_al_(2019)" | QTLID4.LI$Group3=="Current_Study" | QTLID4.LI$Group3=="Clark_et_al_(2016)",]
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)
str(QTLID4.LI)
QTLID5.LI <- QTLID4.LI[,c(1:3,6,10)]
str(QTLID5.LI)
levels(as.factor(QTLID5.LI$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
levels(bands.bc$Group2)
bands.bc1 <- bands.bc[bands.bc$Group2=="Cml"| bands.bc$Group2=="IntL"|bands.bc$Group2=="height"|bands.bc$Group2=="tallest_stem",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc)
bands.bc1 <- droplevels(bands.bc1)
bands.bc2 <- bands.bc1[,c(1,4:5,10,12)]
#bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc2)
levels(as.factor(bands.bc2$Chrom))
####using ID 

chromPlot(gaps=myGM, bands=bands.bc2, stat=QTLID5.LI, statCol="Value",statName="Value")

####using the CmDW
####
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1CmDW.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4.LI <- QTLID4[QTLID4$Group2=="CmDW"| QTLID4$Group2=="CmDW/V", ]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)

#QTLID4.LI <- QTLID4.LI[QTLID4.LI$Group3=="Clark_et_al_(2019)" | QTLID4.LI$Group3=="Current_Study",]
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)
str(QTLID4.LI)
QTLID5.LI <- QTLID4.LI[,c(1:3,6,11)]
str(QTLID5.LI)
levels(as.factor(QTLID5.LI$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
levels(bands.bc$Group2)
bands.bc1 <- bands.bc[bands.bc$Group2=="CmDW"| bands.bc$Group2=="CmDW/V",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc)
bands.bc1 <- droplevels(bands.bc1)
bands.bc2 <- bands.bc1[,c(1,4:5,10,12)]
#bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc2)
levels(as.factor(bands.bc2$Chrom))
####using ID 

chromPlot(gaps=myGM, bands=bands.bc2, stat=QTLID5.LI, statCol="Value",statName="Value")



####using the Yld
####
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1Yld.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4.LI <- QTLID4[QTLID4$Group2=="Yld"| QTLID4$Group2=="SDW"|QTLID4$Group2=="YPP", ]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)

#QTLID4.LI <- QTLID4.LI[QTLID4.LI$Group3=="Clark_et_al_(2019)" | QTLID4.LI$Group3=="Current_Study"| QTLID4.LI$Group3=="Zhao_et_al_(2013)" ,]
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)
str(QTLID4.LI)
QTLID5.LI <- QTLID4.LI[,c(1:3,6,11)]
str(QTLID5.LI)
levels(as.factor(QTLID5.LI$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
levels(bands.bc$Group2)
bands.bc1 <- bands.bc[bands.bc$Group2=="Yld"| bands.bc$Group2=="SDW"|bands.bc$Group2=="YPP",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc)
bands.bc1 <- droplevels(bands.bc1)
bands.bc2 <- bands.bc1[,c(1,4:5,10,12)]
#bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc2)
levels(as.factor(bands.bc2$Chrom))
####using ID 
chromPlot(gaps=myGM, bands=bands.bc2, stat=QTLID5.LI, statCol="Value",statName="Value")




####using the OWA
####
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4.LI <- QTLID4[QTLID4$Group2=="OWA"| QTLID4$Group2=="OWA1"|QTLID4$Group2=="OWA2"| QTLID4$Group2=="OWA3"|
                      QTLID4$Group2=="OWA4"|QTLID4$Group2=="OWA5"|QTLID4$Group2=="OWA6"|QTLID4$Group2=="OWA8"|QTLID4$Group2=="OWA18"|QTLID4$Group2=="OWA19", ]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)

#QTLID4.LI <- QTLID4.LI[QTLID4.LI$Group3=="Clark_et_al_(2019)" | QTLID4.LI$Group3=="Current_Study"| QTLID4.LI$Group3=="Zhao_et_al_(2013)" ,]
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)
str(QTLID4.LI)
QTLID5.LI <- QTLID4.LI[,c(1:3,6,11)]
str(QTLID5.LI)
levels(as.factor(QTLID5.LI$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
levels(bands.bc$Group2)
bands.bc1 <- bands.bc[bands.bc$Group2=="OWA"| bands.bc$Group2=="OWA1"|bands.bc$Group2=="OWA2"|
                        bands.bc$Group2=="OWA3"|bands.bc$Group2=="OWA4"|
                        bands.bc$Group2=="OWA5"|bands.bc$Group2=="OWA6"|bands.bc$Group2=="OWA8"|
                        bands.bc$Group2=="OWA18"|bands.bc$Group2=="OWA19",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc)
bands.bc1 <- droplevels(bands.bc1)
bands.bc2 <- bands.bc1[,c(1,4:5,10,12)]
#bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc2)
levels(as.factor(bands.bc2$Chrom))
####using ID 
chromPlot(gaps=myGM, bands=bands.bc2, stat=QTLID5.LI, statCol="Value",statName="Value")



###flowering time 

####
QTLID4 <- read.csv("Chapter2UP/QTLID1up.org1.csv")
str(QTLID4)
levels(QTLID4$Group2)
levels(QTLID4$Group3)
QTLID4.LI <- QTLID4[QTLID4$Group2=="Flowering"| QTLID4$Group2=="Heading time", ]
#write.csv(QTLID4,file="Chapter2UP/QTLID4BC.csv")
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)

#QTLID4.LI <- QTLID4.LI[QTLID4.LI$Group3=="Clark_et_al_(2019)" | QTLID4.LI$Group3=="Current_Study"| QTLID4.LI$Group3=="Zhao_et_al_(2013)" ,]
QTLID4.LI <- droplevels(QTLID4.LI)
levels(QTLID4.LI$Group2)
levels(QTLID4.LI$Group3)
str(QTLID4.LI)
QTLID5.LI <- QTLID4.LI[,c(1:3,6,11)]
str(QTLID5.LI)
levels(as.factor(QTLID5.LI$Chrom))


###import the bands dataset 
bands.bc <- read.csv("Chapter2UP/Makingbandsgroup1UP.csv")
str(bands.bc)
levels(bands.bc$Group3)
levels(bands.bc$Group2)
bands.bc1 <- bands.bc[bands.bc$Group2=="Flowering"| bands.bc$Group2=="Heading time",]
#bands.bc <- bands[bands$Colors=="Red",]
str(bands.bc)
bands.bc1 <- droplevels(bands.bc1)
bands.bc2 <- bands.bc1[,c(1,4:5,10,12)]
#bands.bc1 <- bands.bc1[order(bands.bc1$Chrom,bands.bc1$Start,bands.bc1$End),]
str(bands.bc2)
levels(as.factor(bands.bc2$Chrom))
####using ID 
chromPlot(gaps=myGM, bands=bands.bc2, stat=QTLID5.LI, statCol="Value",statName="Value")

###this one for groups
QTL <- read.csv("Chapter2UP/ChromPlot.current.study.twogroups5.csv")  
str(QTL)
QTL1 <- QTL[QTL$Group2=="BC"| QTL$Group2=="CCBB"| QTL$Group2=="Bcirc"| QTL$Group2=="CC"| QTL$Group2=="Ccirc",]
QTL3 <- QTL1[,c(1:3,6,4)]
str(QTL3)
levels(QTL$Group2)
#QTL1 <- QTL[QTL$Group2=="BC"| QTL$Group2=="CCBB"| QTL$Group2=="Bcirc",]
chromPlot(gaps=myGM, segment=QTL3, noHist=TRUE, chr=c(1:11),figCols=3, legChrom=11,stack=TRUE)







