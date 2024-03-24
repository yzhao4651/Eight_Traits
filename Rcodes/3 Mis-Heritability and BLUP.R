###import the normalited data，traits1718normalited5 includes
###import the normalited data
normadata <- read.csv("data/traits1718normalited5.csv",na.strings = c("",".","NA"),row.names=1)
###check the data format
str(normadata)
###change the format of the several variables 
normadata$GS <- as.numeric(as.character(normadata$GS))
normadata$SRD <- as.numeric(as.character(normadata$SRD))
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.OWA.R")
Herit <- Heritability(5,38,normadata)
write.csv(Herit, file = "data/heritabilityall5.csv", na = ".") ### this one is for all of the traits for this project. 

###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.R")
ranefvalueall <- ranefvalue(5,38,normadata)
#write.csv(ranefvalueall,file="data/ranefvalueallcheck.csv")

write.csv(ranefvalueall, file = "data/ranefvalueall.csv", row.names = T, na = ".")
