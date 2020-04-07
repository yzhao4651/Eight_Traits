###import the data
###import the data
#normadata <- read.csv("~/Documents/whole traits/traits1718normalited1.csv",na.strings = c("",".","NA"))
#normadata <- read.csv("data/traits1718normalited1.csv",na.strings = c("",".","NA"),row.names=1)
###check the data format
#str(normadata)

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
write.csv(Herit, file = "data/heritabilityall5.csv", na = ".")

###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.R")
ranefvalueall <- ranefvalue(5,38,normadata)
#write.csv(ranefvalueall,file="data/ranefvalueallcheck.csv")

write.csv(ranefvalueall, file = "data/ranefvalueall.csv", row.names = T, na = ".")

###only for OWA 2018,2019
###only for OWA 2018, 2019
#normadata <- read.csv("~/Documents/whole traits/traits1718normalited1.csv",na.strings = c("",".","NA"))
normadata <- read.csv("data/traits1718normalited3.csv",na.strings = c("",".","NA"),row.names=1)
normadata <- normadata[!(normadata$Year=="2017"),]
###check the data format
normadata <- normadata[2:5]
str(normadata)
###change the format of the several variables 
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.R")
Herit <- Heritability(4, 4, normadata)
write.csv(Herit, file = "data/heritabilitOWA1819.csv", na = ".")
###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.OWA.1819.R")
ranefvalueall <- ranefvalue(4, 4, normadata)
write.csv(ranefvalueall, file = "data/ranefvalueOWA1819.csv", na = ".")


#normadata <- read.csv("~/Documents/whole traits/traits1718normalited1.csv",na.strings = c("",".","NA"))
normadata <- read.csv("data/traits1718normalited3.csv",na.strings = c("",".","NA"),row.names=1)
normadata <- normadata[!(normadata$Year=="2019"),]
###check the data format
normadata <- normadata[2:5]
str(normadata)
###change the format of the several variables 
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.R")
Herit <- Heritability(4, 4, normadata)
write.csv(Herit, file = "data/heritabilitOWA1817T.csv", na = ".")
###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.R")
ranefvalueall <- ranefvalue(4, 4, normadata)
write.csv(ranefvalueall, file = "data/ranefvalueOWA1817T.csv", na = ".")





#normadata <- read.csv("~/Documents/whole traits/traits1718normalited1.csv",na.strings = c("",".","NA"))
normadata <- read.csv("data/traits1718normalited3.csv",na.strings = c("",".","NA"),row.names=1)
#normadata <- normadata[!(normadata$Year=="2019"),]
###check the data format
normadata <- normadata[2:5]
str(normadata)
###change the format of the several variables 
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.R")
Herit <- Heritability(4, 4, normadata)
write.csv(Herit, file = "data/heritabilitOWA181719.csv", na = ".")
###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.OWA.R")
ranefvalueall <- ranefvalue(4, 4, normadata)
write.csv(ranefvalueall, file = "data/ranefvalueOWA181719.csv", na = ".")


normadata <- read.csv("data/traits1718normalited3.csv",na.strings = c("",".","NA"),row.names=1)
normadata.17 <- normadata[normadata$Year=="2017",]
normadata.17 <- droplevels(normadata.17)
## OWA will be separated into two years. 
### checking each of year 2018
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.OWA.eachyear.R")
Herit17 <- Heritability(4, 4, normadata.17)
write.csv(Herit17, file = "data/heritabilitOWA17.csv", na = ".")
###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.OWAEachYr.R")
ranefvalueall17 <- ranefvalue(4, 4, normadata.17)
write.csv(ranefvalueall17, file = "data/ranefvalueOWA17.csv", na = ".")

## OWA will be separated into two years. 
### checking each of year 2018
normadata.18 <- normadata[normadata$Year=="2018",]
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.OWA.eachyear.R")
Herit18 <- Heritability(4, 4, normadata.18)
write.csv(Herit18, file = "data/heritabilitOWA18.csv", na = ".")
###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.OWAEachYr.R")
ranefvalueall18 <- ranefvalue(4, 4, normadata.18)
write.csv(ranefvalueall18, file = "data/ranefvalueOWA18.csv", na = ".")


### checking each of year 2019
normadata.19 <- normadata[normadata$Year=="2019",]
write.csv(normadata.19,file="data/OWAdata19.csv")
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.OWA.eachyear.R")
Herit19 <- Heritability(4, 4, normadata.19)
write.csv(Herit19, file = "data/heritabilitOWA19.csv", na = ".")
###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.OWAEachYr.R")
ranefvalueall19 <- ranefvalue(4, 4, normadata.19)
write.csv(ranefvalueall19, file = "data/ranefvalueOWA19.csv", na = ".")
### merge the dataset 18 and 19 together 
OWA17 <- read.csv("data/ranefvalueOWA17.csv")
colnames(OWA17)[colnames(OWA17)=="x"] <- "OWA.17"
OWA18 <- read.csv("data/ranefvalueOWA18.csv")
colnames(OWA18)[colnames(OWA18)=="x"] <- "OWA.18"
OWA19 <- read.csv("data/ranefvalueOWA19.csv")
colnames(OWA19)[colnames(OWA19)=="x"] <- "OWA.19"
OWA1817ranef <- merge(x=OWA17, y=OWA18, by="X")
ranefvalueOWA1817T <- read.csv("data/ranefvalueOWA1817T.csv")
colnames(ranefvalueOWA1817T)[colnames(ranefvalueOWA1817T)=="OWA"] <- "OWA.1718"
OWAall <- merge(x=ranefvalueOWA1817T, y=OWA1817ranef, by="X")
write.csv(OWAall, file="data/OWAallranef.csv")
###import the data of the 


### checking each of year 2019,2018
normadata <- read.csv("data/traits1718normalited3.csv",na.strings = c("",".","NA"),row.names=1)
normadata.1819 <- normadata[normadata$Year=="2019"| normadata$Year=="2018",]
normadata.1819 <- droplevels(normadata.1819)
write.csv(normadata,file="data/OWAdata1819.csv")
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.OWA.eachyear.R")
Herit <- Heritability(4, 4, normadata.1819)
write.csv(Herit, file = "data/heritabilitOWA1819.csv", row.names = T, na = ".")
###get the heritablity of all the traits using the function in the Function


###only for flowering traits
###only for flowering traits
normadata <- read.csv("data/flowpc1.csv",na.strings = c("",".","NA"))
###check the data format
str(normadata)
###change the format of the several variables 
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.R")
Herit <- Heritability(4, 7, normadata)
write.csv(Herit, file = "data/heritabilityflopc1.csv", na = ".")