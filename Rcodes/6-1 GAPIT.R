####GAPIT 
###imported all the data need for GAPIT
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
myYOWA <- read.csv("data/myYimputedSNP19OWA2.csv")
myYSRD <- read.csv("data/myYimputedSNP19SRD.csv")
myY1 <- read.csv("data/myYimputedSNP19.csv")
myY1 <- myY[,-c(2,18)]
myY <- Reduce(function(x,y) merge(x,y,all.x=TRUE),list(myY1,myYOWA,myYSRD))
#myKI <- read.csv("data/kmatrix.df.csv", head = FALSE)
myGD <- read.csv("data/myGDimputedSNP19.csv")
myGM <- read.csv("data/myGMimputedSNP19.csv")
myQ<- read.csv("data/myQimputedSNP19.csv")

###install the packaged need for GAPIT
source("http://www.bioconductor.org/BiocManager.R")
install.packages("BiocManager")
BiocManager::install("multtest")
install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("ape")
install.packages("EMMREML")
install.packages("scatterplot3d") #The downloaded link at: http://cran.r-project.org/package=scatterplot3d
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

###########Using CMLM SUPPER with CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/CMLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/CMLMSUPPERCV0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  SNP.MAF =0.05,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)
###########Using MLM SUPPER with CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/MLMSUPPERCV0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  SNP.MAF =0.05,
  kinship.cluster=c("average"), kinship.group=c("Mean"), 
  sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)

####Enhanced Compression with CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/ECMMLCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/ECMMLCV0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  SNP.MAF =0.05,
  kinship.cluster=c("average"), kinship.group=c("Mean"), 
  group.from=100,
  group.to=36088, 
  group.by=10
)

###########"GLM","MLM","MLMM","FarmCPU" with CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/GLMMLMMLMMCV0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ, 
  SNP.MAF =0.05,
  model=c("GLM","MLM","MLMM")
)




###this without CV
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/CMLMSUPPER0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  SNP.MAF =0.05,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)

setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/MLMSUPPER0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  SNP.MAF =0.05,
  sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)

####Enhanced Compression without CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/ECMMLNOCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/ECMML0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  SNP.MAF =0.05,
  kinship.cluster=c("average"), kinship.group=c("Mean"), 
  group.from=100,
  group.to=36088, 
  group.by=10
)


  
###########"GLM","MLM","MLMM","FarmCPU" without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/GLMMLMMLMM0.05")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/GLMMLMMLMM0.05")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    SNP.MAF =0.05,
    model=c("GLM","MLM","MLMM")
  )
  