####import the data####
####import the data####
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
source("Function/0 Qualdat import.R")
qualdat <- read_qualdat("data/trait1718.3.16.19withoutliers1.csv")
###write out the data set 
write.csv(qualdat,file="~/Documents/whole traits/trait1718SAS1.csv")
write.csv(qualdat,file="data/trait1718SAS1.csv")
####check the data formate
str(qualdat)
### get the chaper2 data sets 
qualdat.chap2 <- qualdat[,c(10:18)]

##get correlation coefficient
qualdat.chap2.coefficient <- cor(qualdat.chap2,use="pairwise")
write.csv(qualdat.chap2.coefficient,file="Charpter2/qualdat.chap2.coefficient.csv")
##import 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
qualdat <- read.csv("data/trait1718SAS1.csv", na.strings = c("",".","NA"),row.names = 1)
str(qualdat)
###seperate dataset using year in order to use later on 
qualdat.17 <- subset(qualdat,qualdat$Year=="2017")
str(qualdat.17)
qualdat.18 <- subset(qualdat,qualdat$Year=="2018")
str(qualdat.18)
#####data summary
#####data summary
####data summary Possible functions used in sapply include mean, sd, var, min, max, median, range, and quantile
install.packages("pastecs")
library(pastecs)
options(scipen=100)
options(digits=4)
stat.desc(qualdat)
summary <-stat.desc(qualdat)
write.csv(summary,file="data/summary.csv",
          eol = "\n", na = "NA")
##for the chapter2 data set 
install.packages("pastecs")
library(pastecs)
options(scipen=100)
options(digits=4)
stat.desc(qualdat.chap2)
summary <-stat.desc(qualdat.chap2)
write.csv(summary,file="data/chapter2.summary.csv",
          eol = "\n", na = "NA")

library(pastecs)
options(scipen=100)
options(digits=4)
stat.desc(qualdat)
summary <-stat.desc(qualdat.17)
write.csv(summary,file="data/summary17.csv",
          eol = "\n", na = "NA")

library(pastecs)
options(scipen=100)
options(digits=4)
stat.desc(qualdat)
summary <-stat.desc(qualdat.18)
write.csv(summary,file="data/summary18.csv",
          eol = "\n", na = "NA")
####scatterplot
source("Function/Scatterplot.R")
###flowering traits with 2017 and 2018
###days as measured method 
str(qualdat)
pdf(paste("Day", 1,".pdf",sep="")) 
Scatterplot(qualdat[,6:9])  

dev.off() 
###flowering traits with 2017 and 2018
###week as measured method 
pdf(paste("Week", 1,".pdf",sep="")) 
Scatterplot(qualdat[,27:30])  
dev.off() 
###flowering traits with 2017 and 2018
###group week as mesured method 
pdf(paste("Group Week", 1,".pdf",sep="")) 
Scatterplot(qualdat[,31:34])  
dev.off() 
###flowering traits with 2017 and 2018
###month as measured method 
source("Function/Scatterplot.month.R")
pdf(paste("Month", 1,".pdf",sep="")) 
Scatterplot.month(qualdat[,35:38])  
dev.off() 


###flowering traits with 2017 and 2018
pdf(paste("Scatterplot of flowering traits", 2,".pdf",sep="")) 
Scatterplot(qualdat[,6:9])  
dev.off() 
###flowering traits with 2017 and 2018
pdf(paste("Scatterplot of flowering traits", 2,".pdf",sep="")) 
Scatterplot(qualdat[,6:9])  
dev.off() 
###flowering traits with 2017
pdf(paste("Scatterplot of 2017 flowering traits", 2,".pdf",sep="")) 
Scatterplot(qualdat.17[,6:9])  
dev.off()
###flowering traits with 2017
pdf(paste("Scatterplot of 2018 flowering traits", 2 ,".pdf",sep="")) 
Scatterplot(qualdat.18[,6:9])  
dev.off()
####other traits with 2017 years
###check the format of the data 

pdf(paste("Scatterplot of 2017 traits", 2 ,".pdf",sep="")) 
Scatterplot(qualdat.17[,c(10:17,19:21,24)])  
dev.off() 
####other traits with 2018 years
###check the format of the data 

pdf(paste("Scatterplot of 2018 traits", 2,".pdf",sep="")) 
Scatterplot(qualdat.18[,c(10:23,25:26)])  
dev.off() 
####other traits with 2017, 2018 years

pdf(paste("Scatterplot of 1718 traits", 2 ,".pdf",sep="")) 
Scatterplot(qualdat[,c(10:17,19:21)])  
dev.off() 
###for charpter2
pdf(paste("Scatterplot of charpter2 traits", 2 ,".pdf",sep="")) 
Scatterplot(qualdat.chap2)  
dev.off()
#####BoxPlot VS Replication
#####BoxPlot VS Replication

source("Function/plotboxFunc.R")
pdf(paste("Boxplot of all traits VS Replication", 2 ,".pdf",sep=""))
plotboxFunc(5,38,qualdat)
dev.off()
#####BoxPlot VS Year
#####BoxPlot VS Year

###select the dataset with both years

qualdat1718 <- qualdat[, c(2,4:17,19:21,27:38)]
source("Function/BoxplotFunVSYear.R")
pdf(paste("Boxplot of all traits VS Year", 2 ,".pdf",sep=" "))
BoxplotFunVSYear(3,30,qualdat1718)
dev.off()
####qq plot and histogram 
####qq plot and histogram 
source("Function/histqq.R")
str(qualdat)
qualdat$OWA <- as.numeric(as.character(qualdat$OWA))
qualdat$CmN. <- as.numeric(as.character(qualdat$CmN.))
qualdat$FNsmall <- as.numeric(as.character(qualdat$FNsmall))
pdf(paste("Q-Q normal plot and Histogram", 2,".pdf",sep=""))
histqq(qualdat[5:38])
dev.off()
###only for flowering traits 
pdf(paste("Q-Q normal plot and Histogram", 2,".pdf",sep=""))
histqq(qualdat[5:38])
dev.off()

###only for flowering traits 
pdf(paste("Q-Q normal plot and Histogram of qualdat.chap2", 1,".pdf",sep=""))
histqq(qualdat.chap2)
dev.off()

###only for flowering traits 
pdf(paste("Q-Q normal plot and Histogram", 1,".pdf",sep=""))
histqq(qualdat[,c(6:9,27:38)])
dev.off()
####Normality
####Normality
#### Boxcox function R codes from online to get lambda

source("Function/bcplot function.txt")
bcplot(na.omit(qualdat[5:38]))
pdf(paste("lambda", 2 ,".pdf",sep=""))####question: only produce one image --> Did we solve this earlier?  seems to product all images now
par(mar=rep(2,4))
par(mfrow=c(4,4))
out_start=5
out_end=38
# set up empty data frame to contain lambda values
lda <- data.frame(row.names = colnames(qualdat)[out_start:out_end],
                  lambda = rep(NA_real_, out_end - out_start + 1))
for (i in out_start:out_end){ 
  trtname <- colnames(qualdat)[i]
  lda[trtname, 1] <- bcplot(na.omit(qualdat[i])) # put result of bcplot into data frame
} 
dev.off()
write.csv(lda,file="data/lda.csv")
###question: how i can get a csv file with column name with lambda?
### Unfortunately the bcplot function does not return any value, so it can't 
### be used programatically to send values to a vector that can be written to a
### file.  You could possible edit the function with a statement like
### `return(lambda.hat)`
###????question below
###I tried to add this lines "return(lambda.hat)", but it does not work, I copied the output into CSV 
###and then import it. I does not take me a lot of time, but it is better to return it as a csv or excel. 
### --> Lindsay's answer
### bcplot was returning lambda.hat, but you didn't update your script to save that value
### anywhere.  I have edited it.

### get the lambdal for the variable need do data transformation 
##quesiton: trying to write a loop but it does not work, do you have any idea?
## The loop was not constructed correctly, since there were not curly brackets
## after it.  I have fixed it, and also indexed by trait name to make things
## easier.
####This one from Lindsay's contribution is GREAT, Thank you so much. 

lda<- read.csv("data/lda.csv",row.name=1)
#lda <- read.csv(file.path(mywd, "lambda1.csv"), row.name=1)

bc1 <- function(x, lda){ # function to transform data after lambda is determined
  for(trait in rownames(lda)){
    while(min(x[[trait]], na.rm = TRUE) <= 0){
      x[[trait]] <- x[[trait]] + 1 # positive numbers required
    }
    if(lda[trait,] == 0){
      x[[trait]] <- log(x[[trait]])
    } else {
      x[[trait]] <- (x[[trait]]^lda[trait,] - 1)/lda[trait,]
    }
  }
  return(x)
}
source("Function/bc1.R")
qualdat_BC <- bc1(qualdat, lda)

####write the data out 
### this one is for owa
write.csv(qualdat_BC, file = "data/traits1718normalited3.csv",row.names = T, na = ".")
## this one is for all of tratis 
write.csv(qualdat_BC, file = "data/traits1718normalited5.csv",row.names = T, na = ".")
#write.csv(qualdat_BC, file = "~/Documents/whole traits/traits1718normalited2.csv", row.names = T, na = ".")


