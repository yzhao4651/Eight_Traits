####GAPIT 
###imported all the data need for GAPIT
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY1 <- read.csv("data/myYimputedSNP19.csv")
myYOWA <- read.csv("data/myYimputedSNP19OWA.csv")
myYSRD <- read.csv("data/myYimputedSNP19SRD.csv")
myY1 <- myY1[,-c(2,18)]
myY <- Reduce(function(x,y) merge(x,y,all.x=TRUE),list(myY1,myYOWA,myYSRD))
str(myY)
names(myY)
myGD <- read.csv("data/myGDimputedSNP19.csv")
myGM <- read.csv("data/myGMimputedSNP19.csv")
myQ<- read.csv("data/myQimputedSNP19.csv")
myQ3<- read.csv("data/myQimputedSNP19.csv", row.names = 1)

###only for FarmCPU
###only for FarmCPU
###only for FarmCPU

  install.packages("bigmemory")
  install.packages("biganalytics")
  library("bigmemory")
  library("biganalytics")
  library("bigmemory")
  library("biganalytics")
  library("compiler") #this library is already installed in R
  source("http://zzlab.net/GAPIT/gapit_functions.txt")
  source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
  str(myY)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/FarmCPUCV0.05")
out_start=2
out_end=38
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/FarmCPU0.05")
out_start=2
out_end=38
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 


###only for FarmCPU
###only for FarmCPU
###only for FarmCPU

install.packages("bigmemory")
install.packages("biganalytics")
library("bigmemory")
library("biganalytics")
library("bigmemory")
library("biganalytics")
library("compiler") #this library is already installed in R
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
str(myY)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYf.106.4185.csv")
str(myY)
myGD <- read.csv("data/myGDf.106.4185.csv")
myGM <- read.csv("data/myGMf.106.4185.csv")
myQ3<- read.csv("data/myQf.106.4185.csv",row.names = 1)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/f106FarmCPUCV0.05")
out_start=2
out_end=24
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/f106FarmCPU0.05")
out_start=2
out_end=24
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYf.116.3098.csv")
myGD <- read.csv("data/myGDf.116.3098.csv")
myGM <- read.csv("data/myGMf.116.3098.csv")
myQ3<- read.csv("data/myQf.116.3098.csv",row.names = 1)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/f116FarmCPUCV0.05UP")
out_start=2
out_end=24
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  
###already run
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/f116FarmCPU0.05UP")
out_start=2
out_end=24
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 


setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYf.96.4814.csv")
myGD <- read.csv("data/myGDf.96.4814.csv")
myGM <- read.csv("data/myGMf.96.4814.csv")
myQ3<- read.csv("data/myQf.96.4814.csv",row.names = 1)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/f96FarmCPUCV0.05")
out_start=2
out_end=24
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/f96FarmCPU0.05")
out_start=2
out_end=24
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 
####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
####becasue i saw the example does not use the PCA, but CV in the manual 
##########Using MLM SUPPER with CV
install.packages("bigmemory")
install.packages("biganalytics")
library("bigmemory")
library("biganalytics")
library("bigmemory")
library("biganalytics")
library("compiler") #this library is already installed in R
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.106.4202.csv")
str(myY)
myGD <- read.csv("data/myGDC.106.4202.csv")
myGM <- read.csv("data/myGMC.106.4202.csv")
myQ3<- read.csv("data/myQC.106.4202.csv",row.names = 1)

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/C106FarmCPUCV0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/C106FarmCPU0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    method.bin="optimum",
    threshold.output=1,
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.116.3293.csv")
myGD <- read.csv("data/myGDC.116.3293.csv")
myGM <- read.csv("data/myGMC.116.3293.csv")
myQ3<- read.csv("data/myQC.116.3293.csv",row.names = 1)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/C116FarmCPUCV0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/C116FarmCPU0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 


setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.124.2560.csv")
str(myY)
myGD <- read.csv("data/myGDC.124.2560.csv")
myGM <- read.csv("data/myGMC.124.2560.csv")
myQ3<- read.csv("data/myQC.124.2560.csv",row.names = 1)

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/C124FarmCPUCV0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/C124FarmCPU0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 
####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
####becasue i saw the example does not use the PCA, but CV in the manual 
##########Using MLM SUPPER with CV

##this one for OWA
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYO.102.4322.csv")
str(myY)
myGD <- read.csv("data/myGDO.102.4322.csv")
str(myGD)
myGM <- read.csv("data/myGMO.102.4322.csv")
str(myGM)
myQ3<- read.csv("data/myQO.102.4322.csv",row.names = 1)
str(myQ3)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/O102FarmCPUCV0.05")
out_start=2
out_end=3
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/O102FarmCPU0.05")
out_start=2
out_end=3
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}
####GAPIT 
###imported all the data need for GAPIT
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYO.112.3450.csv")
str(myY)
myGD <- read.csv("data/myGDO.112.3450.csv")
str(myGD)
myGM <- read.csv("data/myGMO.112.3450.csv")
str(myGM)
myQ3<- read.csv("data/myQO.112.3450.csv",row.names = 1) 
str(myQ)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/O112FarmCPUCV0.05")
out_start=2
out_end=3
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/O112FarmCPU0.05")
out_start=2
out_end=3
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYO.122.2646.csv")
str(myY)
myGD <- read.csv("data/myGDO.122.2646.csv")
str(myGD)
myGM <- read.csv("data/myGMO.122.2646.csv")
str(myGM)
myQ3<- read.csv("data/myQO.122.2646.csv",row.names = 1)
str(myQ3)

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/O122FarmCPUCV0.05")
out_start=2
out_end=3
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/O122FarmCPU0.05")
out_start=2
out_end=3
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}

