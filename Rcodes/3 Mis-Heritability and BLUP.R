###import the normalited data，traits1718normalited_Yld includes all of 38 traits for this study only 8 traits were used for this article
###import the normalited data,traits1718normalited_Yld includes all of 38 traits for this study only 8 traits were used for this article
normadata <- read.csv("data/traits1718normalited_Yld.csv",na.strings = c("",".","NA"),row.names=1)
###check the data format
str(normadata)
###change the format of the several variables 
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)

##change the colnames
colnames(normadata)[which(names(normadata)=="CmDW")] <- "CmDW_g"
colnames(normadata)[which(names(normadata)=="Cml")] <- "Cml_cm"
colnames(normadata)[which(names(normadata)=="CmD_BI")] <- "CmD_BI_mm"
colnames(normadata)[which(names(normadata)=="CmD_LI")] <- "CmD_LI_mm"
colnames(normadata)[which(names(normadata)=="CmN")] <- "CmN."
colnames(normadata)[colnames(normadata)=="Bcirc_cm"] <- "Bcirc"
colnames(normadata)[colnames(normadata)=="CCirc_cm"] <- "CCirc"

###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.R")
ranefvalueall <- ranefvalue(5,17,normadata)

write.csv(ranefvalueall, file = "data/ranefvalueChapter2.csv", row.names = T, na = ".")
