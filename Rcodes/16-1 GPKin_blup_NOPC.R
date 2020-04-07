## testing out genomic prediction for Maria, using kin.blup and not kinship.BlUP

#load("data/160324EMimputedSNP_Msi.RData")
#names(myA.EM.Msi) # $A is the kinship matrix
#singlesiteBLUPs <- read.csv("data/170106BLUPs_to_redo_for_GS_and_GWAS.csv", row.names = 1, header = TRUE)


###import the dataset with phenotype 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYimputedSNP19updated.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDimputedSNP19.csv",row.names=1)
###calculte the K matrix 
# trait to try out
thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GP/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPPC0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:38){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }

    cor(singlesiteBLUPs[theseInd,thisTrait], predValues, use = "complete") # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
}

###import the dataset with phenotype 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYimputedSNP19OWA2.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDimputedSNP19OWA2.csv",row.names=1)
###calculte the K matrix 
# trait to try out
thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GP/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRR0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:2){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues, use = "complete") # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
}


###O.102.4322
###O.102.4322
###this one for nonmissing combinations with 
###import the dataset with phenotype 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYO.102.4322.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDO.102.4322.csv",row.names=1)

###calculte the K matrix 
# trait to try out
thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPOWA102R0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:2){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
}  


####O.112.3450
###O.112.3450
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYO.112.3450.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDO.112.3450.csv",row.names=1)
###calculte the K matrix 
# trait to try out
thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPOWA112R0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:2){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
}

####O.122.2646
###O.122.2646
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYO.122.2646.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDO.122.2646.csv",row.names=1)

###calculte the K matrix 
# trait to try out
thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPOWA122R0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:2){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
}



###f.106.4185
###f.106.4185
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYf.106.4185.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDf.106.4185.csv",row.names=1)
###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPflo106/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:23){
  thisTrait <- traitlist[k,]
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 

###f.116.3098
###f.116.3098
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYf.116.3098.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDf.116.3098.csv",row.names=1)

###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPflo116R0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:23){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 


#####f.96.4814
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYf.96.4814.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDf.96.4814.csv",row.names=1)
###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPflo96R0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:23){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 

####C.106.4202
####C.106.4202
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYC.106.4202.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDC.106.4202.csv",row.names=1)
###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPC106/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRC106R0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:14){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 


####C.116.3293
####C.116.3293
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYC.116.3293.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDC.116.3293.csv",row.names=1)

###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPC116/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRC116R0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:14){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 

#####C.124.2560
#####C.124.2560
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYC.124.2560.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDC.124.2560.csv",row.names=1)
###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPC125/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRC124R0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:14){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 

