
#################################################working for the descriptive statistist for all traits########################
#################################################working for the descriptive statistist for all traits########################

###preparing for importing the data

getwd()##check the path
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")###set up the pathway

###import data and then use the function "0 Qualdat import.R" to organize dataset. 

#The easiest way to get lubridate is to install the whole tidyverse:
install.packages("tidyverse")
#and then get the library "lubridate". intall all of the libraries need. 
install.packages("lubridate")
library(tidyverse)
library(lubridate)

##install the R packages needed later on
install.packages("data.table")
library(data.table)
install.packages("dplyr")
library(dplyr)


source("0 Qualdat import.R") ##using this source R function 

##import the phenotype dataset 
##import the phenotype dataset 
qualdat <- read_qualdat("data/trait1718.3.16.19withoutliers1.csv")



###write out the data set 
write.csv(qualdat,file="~/Documents/whole traits/trait1718SAS1.csv") ##saved a dataset can be used for SAS software for descriptive statistics and saved them under whole traits folder
write.csv(qualdat,file="data/trait1718SAS1.csv") ##saved dataset can be used for SAS software for descriptive analysis 

####check the data formate
str(qualdat)

###R codes for article "Genome-wide association and genomic prediction for biomass traits of Miscanthus sinensis germplasm planted in Alabama"#########
###R codes for article "Genome-wide association and genomic prediction for biomass traits of Miscanthus sinensis germplasm planted in Alabama"#########

####Selected the traits only used for the articles

qualdat.chap2 <- qualdat[,c(1:4,10:18)]
str(qualdat.chap2)

###Change Yld_kg to subsample_Yld_kg (rename it to Yld): if the value > = 1, which will become 1 through the calucatin below:

qualdat.chap2$Yld_kg_sub <- qualdat.chap2$Yld_kg
qualdat.chap2$Yld_kg_sub[qualdat.chap2$Yld_kg_sub >1 ] <- 1

##calculate the Yld
##calculate the Yld
qualdat.chap2$Yld <- (qualdat.chap2$Yld_kg*qualdat.chap2$SDW_kg/qualdat.chap2$Yld_kg_sub)*1000

###write out the chapter datasets 

write.csv(qualdat.chap2,file="Chapter2UP/qualdat.chap2original.csv")

##get correlation coefficient
##get correlation coefficient
qualdat.chap2.coefficient <- cor(qualdat.chap2,use="pairwise")
write.csv(qualdat.chap2.coefficient,file="Charpter2/qualdat.chap2.coefficient.csv")


###data summary to get the mean, max, min, range, variance, and standard error for Table 4 A 
###data summary to get the mean, max, min, range, variance, and standard error for Table 4 A 

install.packages("pastecs")
library(pastecs)
options(scipen=100)
options(digits=4)
stat.desc(qualdat.chap2)
summary <-stat.desc(qualdat.chap2)
write.csv(summary,file="data/chapter2.summary.csv",
          eol = "\n", na = "NA")

####get the normalized datasets for traits only for 8 triats articles 
####get the normalized datasets for traits only for 8 triats articles 

####get lambda for Normality
####get lambda for Normality

####using bcplot function (Boxcox function R codes from online) to get lambda
####using bcplot function (Boxcox function R codes from online) to get lambda

source("Function/bcplot function.txt")

#bcplot(na.omit(qualdat[5:38]))
##change the columne names
colnames(qualdat.chap2)[which(names(qualdat.chap2)=="CmDW_g")] <- "CmDW"
colnames(qualdat.chap2)[which(names(qualdat.chap2)=="Cml_cm")] <- "Cml"
colnames(qualdat.chap2)[which(names(qualdat.chap2)=="CmD_BI_mm")] <- "CmD_BI"
colnames(qualdat.chap2)[which(names(qualdat.chap2)=="CmD_LI_mm")] <- "CmD_LI"
colnames(qualdat.chap2)[which(names(qualdat.chap2)=="CmN.")] <- "CmN"
colnames(qualdat.chap2)[colnames(qualdat.chap2)=="Bcirc_cm"] <- "Bcirc"
colnames(qualdat.chap2)[colnames(qualdat.chap2)=="CCirc_cm"] <- "CCirc"
#colnames[colnames(qualdat.chap23)=="CmDW_g"] <- "CmDW"

qualdat <- qualdat.chap2[,-c(11,14)]
bcplot(na.omit(qualdat))
pdf(paste("lambda", 2 ,".pdf",sep=""))####question: only produce one image --> Did we solve this earlier?  seems to product all images now
par(mar=rep(2,4))
par(mfrow=c(4,4))
out_start=5
out_end=17
# set up empty data frame to contain lambda values
lda <- data.frame(row.names = colnames(qualdat)[out_start:out_end],
                  lambda = rep(NA_real_, out_end - out_start + 1))
for (i in out_start:out_end){ 
  trtname <- colnames(qualdat)[i]
  lda[trtname, 1] <- bcplot(na.omit(qualdat[i])) # put result of bcplot into data frame
} 
dev.off()

write.csv(lda,file="data/lda.makeup_kg.csv")


lda <- read.csv("data/lda.makeup_kg.csv",row.name=1)
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
qualdat_BC_kg <- bc1(qualdat, lda)
str(qualdat_BC_kg)

####write the data out 
write.csv(qualdat_BC_kg, file = "data/traits1718normalited_Yld.csv", row.names = T, na = ".")


###import the phonetype datasets after data transformation (follow normal distribution)

qualdat.chap2.normal <- read.csv("data/traits1718normalited_Yld.csv")

str(qualdat.chap2.normal)

qualdat.chap2.normal_1 <- qualdat.chap2.normal[,c(2:13,16)]

###import the table with entry and accessions 

Altable <- read.csv("data/Indi.Entry.Access.csv", header = TRUE, na="")
str(Altable)

###import the datasets after data transformed 

## Merging both datasets

#install.packages("plyr")

library("plyr")
all_accession_all <- plyr::join_all(list(Altable,qualdat.chap2.normal_1), by='Entry')
str(all_accession_all)
levels(as.factor(all_accession_all$Accession))

#check_3 <- subset(all_accession_all,all_accession_all$Accession=="PMS-300")

levels(as.factor(all_accession_all$plot.1))


#Remove the number below from the column plot.1 for duplicated or multiple records. 
#"104.1","217.1""222.1" "228.1" "288.1"  "464.1""503.3""514.1"  "528.1" "616.1""639.1""643.1"

library(dplyr)

all_accession_all <- all_accession_all %>%
  filter_at(vars("plot.1"), all_vars(! . %in% c("69.1","104.1","217.1","222.1","228.1","288.1" ,"464.1","503.3","514.1","528.1","616.1","639.1","643.1")))


levels(as.factor(all_accession_all$plot.1))

### Maually change the PMS-14 to PMS-014 and PMS-75 to PMS-075

all_accession_all$Accession <- as.factor(all_accession_all$Accession)
levels(all_accession_all$Accession)[levels(all_accession_all$Accession)=="PMS-14"] <- "PMS-014"
levels(all_accession_all$Accession)[levels(all_accession_all$Accession)=="PMS-75"] <- "PMS-075"
levels(all_accession_all$Accession)[levels(all_accession_all$Accession)=="PMS-9"] <- "PMS-009"
levels(all_accession_all$Accession)[levels(all_accession_all$Accession)=="PMS-38"] <- "PMS-038"

levels(all_accession_all$Accession)


###Import all of the dataset from the Clark.et al 2014 to merge with the data in this study to get their orginal genotypic group
###Genetic structure of Miscanthus sinensis and Miscanthus sacchariflorus in Japan indicates a gradient of bidirectional but asymmetric introgression

##import the dataset from Clark.et al 2014
Orig_1 <- read.csv("data/gcbb12606-sup-0001-datas1_LVC.csv")
str(Orig_1)

levels(as.factor(Orig_1$Genotype))
levels(as.factor(Orig_1$DAPC.group))

str(Orig_1)
names(Orig_1)

Orig_1[4:13] <- round(Orig_1[4:13],1)

names(Orig_1)[1] <- "Individual"

names(Orig_1)

Orig_1$Individual <- as.factor(Orig_1$Individual)

levels(Orig_1$Individual)

levels(Orig_1$Individual)[levels(Orig_1$Individual)=="DC-2010-001-A"] <- "DC-2010-001 A"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="DC-2010-001-E"] <- "DC-2010-001 E"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="DC-2010-001-D"] <- "DC-2009-001 D" 

levels(Orig_1$Individual)[levels(Orig_1$Individual)=="KY-2009-001-C"] <- "KY-2009-001 C"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="KY-2009-001-D"] <- "KY-2009-001 D"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="KY-2009-001-B20-a"] <- "KY-2009-001-B(20)-a"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="KY-2009-001-B20-b"] <- "KY-2009-001-B(20)-b" 
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="KY-2009-001-B20-e"] <- "KY-2009-001-B(20)-e"

levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PA-2010-002-E"] <- "PA-2010-002 E" 
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PA-2010-001-B63-a"] <- "PA-2010-001-B(63)-a" 
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PA-2010-001-B63-c"] <- "PA-2010-001-B(63)-c" 
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PA-2010-002-B3-b"] <- "PA-2010-002-B(3)-b" 
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PA-2010-002-B3-c"] <- "PA-2010-002-B(3)-c" 
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PA-2010-002-B3-d"] <- "PA-2010-002-B(3)-d" 

levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PA-2010-003-B20-b"] <- "PA-2010-003-B(20)-b" 
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PA-2010-003-B20-c"] <- "PA-2010-003-B(20)-c" 


levels(Orig_1$Individual)[levels(Orig_1$Individual)=="VA-2010-001-B40-e"] <- "VA-2010-001-B(40)-e"


levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-001-5-C"] <- "NC-2010-001.5 C"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-001-5-B18-a"] <- "NC-2010-001.5-B(18)-a"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-001-5-B18-b"] <- "NC-2010-001.5-B(18)-b"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-004-E"] <- "NC-2010-004 E" 
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-001-B"] <- "NC-2010-001 B"

levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-003-B50-a"] <- "NC-2010-003-B(50)-a"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-003-B50-b"] <- "NC-2010-003-B(50)-b"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-004-B37-b"] <- "NC-2010-004-B(37)-b"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-004-B37-d"] <- "NC-2010-004-B(37)-d"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-003-D"] <- "NC-2010-003 D"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-003-E"] <- "NC-2010-003 E"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="NC-2010-001-B44-e"] <- "NC-2010-001-B(44)-e"

levels(Orig_1$Individual)[levels(Orig_1$Individual)=="OH-2009-001-B20-b"] <- "OH-2009-001-B(20)-b"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="OH-2009-001-B20-d"] <- "OH-2009-001-B(20)-d"

levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PI230189-US56-0022-03"] <- "PI230189"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PI294602-US64-0004-02"] <- "PI294602"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PI294605-US64-0007-01"] <- "PI294605"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PI230189-US56-0022-03"] <- "PI230189"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PI295764-Gracillimus"] <- "PI295764"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PI302423-US64-0016-03"] <- "PI302423"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PI417947-NG77-022"] <- "PI417947"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="PI423566-M75-062"] <- "PI423566"
levels(Orig_1$Individual)[levels(Orig_1$Individual)=="CANE9233-US47-0011"] <- "CANE9233"

levels(as.factor(Orig_1$DAPC.group))

###import orginal genotypic group using another source from LVC
###import orginal genotypic group using another source from LVC

Orig_2 <- read.csv("data/jexbot135939_file002_LVC.csv")
str(Orig_2)

Orig_3 <- Orig_2[,c(2:3,5:6,11:19)]

str(Orig_3)

Orig_3[3:12] <- round(Orig_3[3:12],1)
str(Orig_3)

levels(as.factor(Orig_3$Individual))

levels(as.factor(Orig_3$DAPC.group))

##remove the DAPC.group if missing 

Orig_3$DAPC.group <- as.factor(Orig_3$DAPC.group)

levels(Orig_3$DAPC.group)[levels(Orig_3$DAPC.group)=="1"] <- "Korea/N China"
levels(Orig_3$DAPC.group)[levels(Orig_3$DAPC.group)=="2"] <- "Yangtze-Qinling"
levels(Orig_3$DAPC.group)[levels(Orig_3$DAPC.group)=="5"] <- "SE China/tropical"
levels(Orig_3$DAPC.group)[levels(Orig_3$DAPC.group)=="7"] <- "Msa"
levels(Orig_3$DAPC.group)[levels(Orig_3$DAPC.group)=="8"] <- "Sichuan"

levels(Orig_3$DAPC.group)[levels(Orig_3$DAPC.group)=="3"] <- "Central Japan"
levels(Orig_3$DAPC.group)[levels(Orig_3$DAPC.group)=="4"] <- "N Japan"
levels(Orig_3$DAPC.group)[levels(Orig_3$DAPC.group)=="6"] <- "S Japan"

levels(as.factor(Orig_3$DAPC.group))



#Index <- which(Orig_3$DAPC.group==1 | Orig_3$DAPC.group==2 | Orig_3$DAPC.group==3 | Orig_3$DAPC.group==4
#|Orig_3$DAPC.group==5  |Orig_3$DAPC.group==6 |Orig_3$DAPC.group==7 |Orig_3$DAPC.group==8)

#Orig_4 <- Orig_3[Index,]
#names(Orig_4)

str(Orig_3)

Orig_4_1 <- Orig_3[,c(1,13,3:12)]

str(Orig_4_1)

Orig_4_1[3:12] <- round(Orig_4_1[3:12],1)

colnames(Orig_4_1)[5] <- "Korea.N.China"

levels(as.factor(Orig_4_1$Individual))

Orig_4_2 <-  Orig_4_1 %>% na.omit(Orig_4_1$DAPC.group)


Orig_4_2$Individual <- as.factor(Orig_4_2$Individual)

levels(Orig_4_2$Individual)

levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="DC-2010-001-A"] <- "DC-2010-001 A"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="DC-2010-001-E"] <- "DC-2010-001 E"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="DC-2010-001-D"] <- "DC-2009-001 D" 

levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="KY-2009-001-C"] <- "KY-2009-001 C"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="KY-2009-001-D"] <- "KY-2009-001 D"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="KY-2009-001-B20-a"] <- "KY-2009-001-B(20)-a"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="KY-2009-001-B20-b"] <- "KY-2009-001-B(20)-b" 
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="KY-2009-001-B20-e"] <- "KY-2009-001-B(20)-e"

levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PA-2010-002-E"] <- "PA-2010-002 E" 
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PA-2010-001-B63-a"] <- "PA-2010-001-B(63)-a" 
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PA-2010-001-B63-c"] <- "PA-2010-001-B(63)-c" 
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PA-2010-002-B3-b"] <- "PA-2010-002-B(3)-b" 
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PA-2010-002-B3-c"] <- "PA-2010-002-B(3)-c" 
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PA-2010-002-B3-d"] <- "PA-2010-002-B(3)-d" 

levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PA-2010-003-B20-b"] <- "PA-2010-003-B(20)-b" 
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PA-2010-003-B20-c"] <- "PA-2010-003-B(20)-c" 

levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="VA-2010-001-B40-e"] <- "VA-2010-001-B(40)-e"

levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-001-5-C"] <- "NC-2010-001.5 C"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-001-5-B18-a"] <- "NC-2010-001.5-B(18)-a"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-001-5-B18-b"] <- "NC-2010-001.5-B(18)-b"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-004-E"] <- "NC-2010-004 E" 
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-001-B"] <- "NC-2010-001 B"

levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-003-B50-a"] <- "NC-2010-003-B(50)-a"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-003-B50-b"] <- "NC-2010-003-B(50)-b"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-004-B37-b"] <- "NC-2010-004-B(37)-b"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-004-B37-d"] <- "NC-2010-004-B(37)-d"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-003-D"] <- "NC-2010-003 D"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-003-E"] <- "NC-2010-003 E"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="NC-2010-001-B44-e"] <- "NC-2010-001-B(44)-e"

levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="OH-2009-001-B20-b"] <- "OH-2009-001-B(20)-b"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="OH-2009-001-B20-d"] <- "OH-2009-001-B(20)-d"

levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PI230189-US56-0022-03"] <- "PI230189"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PI294602-US64-0004-02"] <- "PI294602"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PI294605-US64-0007-01"] <- "PI294605"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PI230189-US56-0022-03"] <- "PI230189"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PI295764-Gracillimus"] <- "PI295764"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PI302423-US64-0016-03"] <- "PI302423"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PI417947-NG77-022"] <- "PI417947"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="PI423566-M75-062"] <- "PI423566"
levels(Orig_4_2$Individual)[levels(Orig_4_2$Individual)=="CANE9233-US47-0011"] <- "CANE9233"

levels(Orig_4_2$Individual)

###both sources create some duplicated genotype.The duolicationsgenotype were produced by manually
###remove some of individula with the same genotypic group using the dataset called duolicationsgenotype

Remove <- read.csv("data/duplicationsgenotype.csv")
str(Remove)
Remove$Var1 <- as.factor(Remove$Var1)

##levels(Remove$Var1)[levels(Remove$Var1)=="PI417947"] <- "PI417947-NG77-022"

Orig_4_3 <- Orig_4_2[!(Orig_4_2$Individual %in% Remove$Var1),]

colnames(Orig_4_3) %in% colnames(Orig_1)

colnames(Orig_1) %in% colnames(Orig_4_3)

##get the dataset by combine both datasets. 

Orig <- unique(rbind(Orig_4_3, Orig_1[,c(1,3:13)]))

levels(as.factor(Orig$DAPC.group))

Orig$DAPC.group <- as.factor(Orig$DAPC.group)


##need to change several indiviuals during the writing format difference

Orig$Individual <- as.factor(Orig$Individual)

levels(Orig$Individual)

##write out the dataset to check the datasets manually 

write.csv(Orig, file="data/check.csv")

##get the unique dataset 
Orig <- unique(Orig)

write.csv(Orig,file="data/orginal_all_LVC_updated.csv") ##write out the original genotype for 


Orig$Individual <-  as.factor(Orig$Individual)

levels(Orig$Individual)

levels(Orig$Individual)

###get the matched genetic group using merge function with cultivars from this study

levels(all_accession_all$Accession) 

levels(Orig$Individual)


Acces.orig_all <- merge(Orig,all_accession_all,by.x=c("Individual"),by.y=c("Accession"),all.y=TRUE)

#Acces.orig <- Orig[Orig$Individual %in% all_accession_all$Accession,] this one also can be used for checking if they match. 
#Acces.orig <- droplevels(Acces.orig)

str(Acces.orig_all)

levels(as.factor(Acces.orig_all$Entry))

##Romved the duplicated cases

Acces.orig_all_sub <- Acces.orig_all[complete.cases(Acces.orig_all[,c("Entry")]),]


factor(Acces.orig_all_sub$DAPC.group)

levels(as.factor(Acces.orig_all_sub$DAPC.group))

write.csv(Acces.orig_all_sub, file="Chapter2UP/Acces.orig_all_sub_check.csv") ##write out the dataset for checking

###add the "Ornamental" for UI12-00005.
#install.packages("tidyr")

library(tidyr)

levels(as.factor(Acces.orig_all_sub$DAPC.group))

Acces.orig_all_sub$DAPC.group  <-  Acces.orig_all_sub$DAPC.group %>% replace_na('Ornamental')

levels(Acces.orig_all_sub$DAPC.group)

all(is.na(Acces.orig_all_sub$DAPC.group))

str(Acces.orig_all_sub)

Acces.orig_all_sub <- droplevels(Acces.orig_all_sub)

Acces.orig_all_sub_unique <- unique(Acces.orig_all_sub)


###check some duplicated values in file datasets and select the freq > 8 to remove 

T <- as.data.frame(table(Acces.orig_all_sub_unique$Individual))
str(T)

T_sub <- subset(T, Freq==16) ###find if that is O and then write out the dataset Acces.orig_all_sub_unique


#write.csv(T_sub,file="data/duplicationsgenotype.csv",row.names = F)

##export the datasets with orgianl genetic type and then double check

write.csv(Acces.orig_all_sub_unique, file="data/Accession.Source.traits_updated_2.csv",row.names = F)

##take only first four columns form 

str(Acces.orig_all_sub_unique)

levels(as.factor(Acces.orig_all_sub$plot.1))

Access.Source.traits_1 <- Acces.orig_all_sub[,c(1:15)]
Access.Source.traits_unique <- unique(Access.Source.traits_1)

str(Access.Source.traits_unique)
nrow(Access.Source.traits_unique)
levels(as.factor(Access.Source.traits_unique$Individual))

write.csv(Access.Source.traits_unique,file="data/Access.Source.traits_unique_up.csv")


### and then remove the M. sacchariflorus and M. ×giganteus ‘Illinois’ 

levels(as.factor(Access.Source.traits_unique$group))

#install.packages("dplyr")
library("dplyr")

Access.Source.traits_unique <- Access.Source.traits_unique %>%
  filter_at(vars("group"), all_vars(! . %in% c("M. ×giganteus","M. sacchariflorus","M. floridulus")))

Access.Source.traits_unique <- droplevels(Access.Source.traits_unique)

levels(as.factor(Access.Source.traits_unique$group))

Access.Source.traits_unique$DAPC.group <- as.factor(Access.Source.traits_unique$DAPC.group)

#levels(Access.Source.traits_unique$DAPC.group)[levels(Access.Source.traits_unique$DAPC.group) == "M. xgiganteus 2x"] <- "M. sinensis var. purpurascens"

levels(Access.Source.traits_unique$DAPC.group)

## produce the genetic group for the section of 2.1 | Plant materials, field trials and phenotyping

table(Access.Source.traits_unique$DAPC.group) 

#Table results for getting their orginal genotypic group for the number in page 6
#Table results for getting their orginal genotypic group for the number in page 6

#Korea/N China   Yangtze-Qinling   SE China/tropical  Sichuan    US naturalized    M. xgiganteus 2x        Ornamental 
#2                19                  15                 7                38                 2                77 


####This is for Figure 2 Scatter plots
####This is for Figure 2 Scatter plots

qualdat_check<- read.csv("data/traits1718normalited_Yld.csv")
str(qualdat_check)
qualdat.chap23 <- qualdat_check[, c(6:11,13,16)]
str(qualdat.chap23)
names(qualdat.chap23)

###change column names
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="CmDW_g")] <- "CmDW"
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="Cml_cm")] <- "Cml"
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="CmD_BI_mm")] <- "CmD_BI"
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="CmD_LI_mm")] <- "CmD_LI"
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="CmN.")] <- "CmN"
colnames(qualdat.chap23)[colnames(qualdat.chap23)=="Bcirc_cm"] <- "Bcirc"
colnames(qualdat.chap23)[colnames(qualdat.chap23)=="CCirc_cm"] <- "CCirc"
names(qualdat.chap23)
qualdat.chap23$CmDW<- as.numeric(qualdat.chap23$CmDW)
qualdat.chap23$Cml<- as.numeric(qualdat.chap23$Cml)
qualdat.chap23$CmD_BI <- as.numeric(qualdat.chap23$CmD_BI)
qualdat.chap23$CmD_LI <- as.numeric(qualdat.chap23$CmD_LI)
qualdat.chap23$CmN <- as.numeric(qualdat.chap23$CmN)
qualdat.chap23$Bcirc <- as.numeric(qualdat.chap23$Bcirc)
qualdat.chap23$CCirc <- as.numeric(qualdat.chap23$CCirc)
qualdat.chap23$Yld <- as.numeric(qualdat.chap23$Yld)
str(qualdat.chap23)

#pdf(paste("Scatterplot of charpter2 usingnormalizedtraits", 2 ,".pdf",sep="")) 
#Scatterplot(qualdat.chap23)  
dev.off()
#install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
pdf(paste("Scatterplot of charpter2 usingnormalizedtraits", 5,".pdf",sep="")) ###saved it as PDF as figure 2
str(qualdat.chap23)
#qualdat.chap23 <- qualdat.chap2[, -c(7:8,10)]
###this one used for scatter plot with 8 traits with normalized dataset
#qualdat.chap23 <- qualdat[, c(6:11,13,16)]
#str(qualdat.chap23)
qualdat.chap23$CmDW<- as.numeric(qualdat.chap23$CmDW)
qualdat.chap23$Cml<- as.numeric(qualdat.chap23$Cml)
qualdat.chap23$CmD_BI <- as.numeric(qualdat.chap23$CmD_BI)
qualdat.chap23$CmD_LI <- as.numeric(qualdat.chap23$CmD_LI)
qualdat.chap23$CmN <- as.numeric(qualdat.chap23$CmN)
qualdat.chap23$Bcirc <- as.numeric(qualdat.chap23$Bcirc)
qualdat.chap23$CCirc <- as.numeric(qualdat.chap23$CCirc)
qualdat.chap23$Yld <- as.numeric(qualdat.chap23$Yld)
colnames(qualdat.chap23)
###change column names
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="CmDW_g")] <- "CmDW"
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="Cml_cm")] <- "Cml"
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="CmD_BI_mm")] <- "CmD_BI"
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="CmD_LI_mm")] <- "CmD_LI"
colnames(qualdat.chap23)[which(names(qualdat.chap23)=="CmN.")] <- "CmN"
colnames(qualdat.chap23)[colnames(qualdat.chap23)=="Bcirc_cm"] <- "Bcirc"
colnames(qualdat.chap23)[colnames(qualdat.chap23)=="CCirc_cm"] <- "CCirc"
#colnames[colnames(qualdat.chap23)=="CmDW_g"] <- "CmDW"
str(qualdat.chap23)
#Scatterplot(qualdat.chap23)
###chart.Correlation 
chart.Correlation(qualdat.chap23, histogram=TRUE,qqplot=TRUE, pch=19, sig.level = 0.00001)
#text(main="The scatter plot, histogram, and correlation coefficient of traits",cex.lab=2)
#dev.off()
#mtext("The scatter plot, histogram, and correlation coefficient of traits", side=3, line=3, cex.lab=6,font=TRUE)
#mtext("CmDW   Cml   CmD_BI  CmD_LI  CmN   Bcirc  Ccirc  Yld ", side=1, line=1, outer=TRUE)
#xlab=("CmDW   Cml   CmD_BI  CmD_LI  CmN   Bcirc  Ccirc  Yld ")
mtext("             CmDW       Cml      CmD_BI    CmD_LI    CmN       Bcirc        Ccirc         Yld ", side=1, line=3, cex.lab=4,padj = 1,adj = 0)
dev.off()

###it also can use this function
Scatterplot <- function(mydata){
  panel.cor <- function(x, y, digits=4, prefix="", cex.cor) 
    
  {
    usr <- par("usr"); on.exit(par(usr)) 
    par(usr = c(0, 1, 0, 1)) 
    r <- abs(cor(x, y, use="pairwise")) 
    txt <- format(c(r, 0.123456789), digits=digits)[1] 
    txt <- paste(prefix, txt, sep="") 
    if(missing(cex.cor)) cex <- 0.9/strwidth(txt) 
    
    test <- cor.test(x, y, use="pairwise")
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.00001, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", " *", ".", " ")) 
    #Signif <- symnum(r, corr = FALSE, na = FALSE, 
    #cutpoints = c(0.3, 0.6, 0.7, 0.8, 0.95),
    #symbols = c("", ".", "*", "**", "***")) 
    
    #symnum(r, cutpoints = c(0.3, 0.6, 0.8, 0.9, 0.95),
    #symbols = c(" ", ".", ",", "+", "*", "B"),
    #abbr.colnames = TRUE)
    text(0.5, 0.5, txt, cex = cex * r) 
    text(.5, .8, Signif, cex = cex, col=2)
  }
  pairs(mydata, lower.panel=panel.smooth, upper.panel=panel.cor)
}


###this one is to get the lsmean for each of accession for Table 4 B 
###this one is to get the lsmean for each of accession  for Table 4 B 

#install.packages(c("ggpubr", "broom", "AICcmodavg"))
library(ggpubr)
library(broom)
library(AICcmodavg)

library(nlme)
library(lme4)
#install.packages("lmerTest")
library(lmerTest)

###bivariates analysis 
###get aova from regression model

str(Acces.orig_all_sub_unique)

levels(as.factor(Acces.orig_all_sub_unique$DAPC.group))

all_accession_1 <- Acces.orig_all_sub_unique[-25]

names(all_accession_1)

levels(as.factor(all_accession_1$DAPC.group))


###using normalized data to do Anova analysis 
###using normalized data to do Anova analysis 

###using the dataset

str(all_accession_1)

all_accession_1$Entry1 <- as.factor(all_accession_1$Entry1)

all_accession_1$DAPC.group <- as.factor(all_accession_1$DAPC.group)

all_accession_1$group <- as.factor(all_accession_1$group)
all_accession_1$Rep <- as.factor(all_accession_1$Rep)
all_accession_1$Year <- as.factor(all_accession_1$Year)

all_accession_1$CmDW <- as.numeric(all_accession_1$CmDW)
all_accession_1$Cml <- as.numeric(all_accession_1$Cml)
all_accession_1$CmD_BI <- as.numeric(all_accession_1$CmD_BI)
all_accession_1$CmD_LI <- as.numeric(all_accession_1$CmD_LI)
all_accession_1$CmN <- as.numeric(all_accession_1$CmN)
all_accession_1$Bcirc <- as.numeric(all_accession_1$Bcirc)
all_accession_1$CCirc <- as.numeric(all_accession_1$CCirc)
all_accession_1$Yld <- as.numeric(all_accession_1$Yld)

str(all_accession_1)

levels(as.factor(all_accession_1$group))
levels(as.factor(all_accession_1$Entry1))
levels(all_accession_1$DAPC.group)

####get the Lsmean 
####get the Lsmean 

library("lme4")
#install.packages("lme4")
str(all_accession_1)
names(all_accession_1)[19]

#install.packages("mvtnorm")
#install.packages("emmeans")
#install.packages("lsmeans")
#install.packages("devtools")
library(devtools)
#remotes::install_github("rvlenth/emmeans")

#install.packages("pbkrtest")
#install.packages("lmerTest")

#all_accession_1$DAPC.group <- factor(all_accession_1$DAPC.group, levels = c("M. ×giganteus","M. floridulus","M. sacchariflorus","M. sinensis","M. sinensis var. condensatus",
#"M. sinensis var. purpurascens","M. sp","M. transmorrisonensis"))
str(all_accession_1)
names(all_accession_1)
levels(all_accession_1$DAPC.group)
#all_accession_1 <- subset(all_accession_1, all_accession_1$DAPC.group!="M. sacchariflorus 2x")

#all_accession_2 <- subset(all_accession_1, all_accession_1$DAPC.group!="M. xgiganteus 3x")

#all_accession_2$DAPC.group <- droplevels(all_accession_2$DAPC.group)

levels(all_accession_1$DAPC.group)

library(lme4)
library("lsmeans")
mz_3 <- list()
mz_4 <- list()
str(all_accession_1)
start=19
end=26
for(i in start:end){
  if(sum(is.na(all_accession_1[,i])) >= 900){
    #Misvarcomp<- lmer(all_accession_1[,i] ~ all_accession_1$group + (1|Rep))
    LM_5 <- lmer(all_accession_1[,i] ~ Entry1 + (1|Rep), data=all_accession_1,na.action=na.exclude)
    #anova(LM_2)
    #drop1(LM_5,test="Chisq")
    #fixef(LM_3)
    #names(summary(LM_3))$coefficients
    #summary(LM_4)$coefficients
    t1 <- data.frame(lsmeans(LM_5, c("Entry1")))
    mz_3[[i]] = data.frame(t1, estimate=rownames(t1), variables= rep(names(all_accession_1[i]), dim(t1)[1]))
  } else {
    #Misvarcomp<- lmer(all_accession_1[,i] ~ (1|Entry) + (1|Rep) + (1|Year))
    LM_6 <- lmer(all_accession_1[,i] ~ Entry1 +(1|Rep)+(1|Year), data=all_accession_1,na.action=na.exclude) ### ANOVA
    t2 <- data.frame(lsmeans(LM_6, c("Entry1")))
    mz_4[[i]] = data.frame(t2, estimate=rownames(t2), variables= rep(names(all_accession_1[i]), dim(t2)[1]))
  }
}
smry_chap2_2 <- plyr::ldply(mz_3, data.frame)
smry_chap2_1 <- plyr::ldply(mz_4, data.frame)
smrylsmean_chap2 <- rbind(smry_chap2_2,smry_chap2_1)

head(smrylsmean_chap2)

levels(as.factor(smrylsmean_chap2$variables))

write.csv(smrylsmean_chap2, file="Chapter2UP/Eighttrait_lsmean.csv",row.names = F)


#### back-transform BLUPs to have real units
#### back-transform BLUPs to have real units
#library(lme4)
#library(arm)
library(lsmeans)
#load("170206phenotype_analysis.RData")

# function for reverse Box-Cox transformation

lda <- read.csv("data/lda.makeup_kg.csv")
lda_1 <-lda[-c(7,9,10,12,13),] 
##change the first column in to row.name

#rownames(lda_1) <- lda_1[,1]
#lda_1[,1] <- NULL
#lda_1

##add lda_1 to the datasets of smrylsmean

str(smrylsmean_chap2)

lsmeanlda <- merge(smrylsmean_chap2, lda_1, by.x=c("variables"),by.y=c("X"), alll.x=TRUE)

str(lsmeanlda)

bcBack <- function(x1, lda){
  ifelse(lda==0, exp(x1), (lda * x1 + 1)^(1/lda) )
}

# back-transform

lsmeanlda$BackLSmeans <- bcBack(lsmeanlda$lsmean, lsmeanlda$lambda)
lsmeanlda$BackLSmeansSE <- bcBack(lsmeanlda$SE, lsmeanlda$lambda)
str(lsmeanlda)

###long to wide 

names(lsmeanlda)

lsmeanlda_1 <- lsmeanlda[,c(1,2,10)]
names(lsmeanlda_1)
str(lsmeanlda_1)

lsmeanlda_2 <- lsmeanlda[,c(1,2,11)]
names(lsmeanlda_2)

library(dplyr)
library(reshape2)
#install.packages("reshape2")
#install.packages("dplyr")

lsmeanlda_Ls <- dcast(lsmeanlda_1, Entry1~ variables, value.var=c("BackLSmeans"))
lsmeanlda_SE<- dcast(lsmeanlda_2, Entry1~ variables, value.var=c("BackLSmeansSE"))

str(lsmeanlda_SE)
str(lsmeanlda_Ls)

names(lsmeanlda_SE)[2] <- "Bcirc_SE"
names(lsmeanlda_SE)[3] <- "CCirc_SE"
names(lsmeanlda_SE)[4] <- "CmD_BI_SE"
names(lsmeanlda_SE)[5] <- "CmD_LI_SE"
names(lsmeanlda_SE)[6] <- "CmDW_SE"
names(lsmeanlda_SE)[7] <- "Cml_SE"
names(lsmeanlda_SE)[8] <- "CmN_SE"
names(lsmeanlda_SE)[9] <- "Yld_SE"

##merge both datasets
Lsmean_se <- merge(lsmeanlda_Ls,lsmeanlda_SE, by=c("Entry1"), all.x=TRUE, all.y=TRUE)

str(Lsmean_se)

##merge the genotipic group 
str(all_accession_1)
names(all_accession_1)
all_accession_2 <- all_accession_1[,c(1:2,13:15)]
Lsmean_genetic <- merge(all_accession_2, Lsmean_se, by=c("Entry1"), all.x=TRUE, all.y=TRUE)
str(Lsmean_genetic)

names(Lsmean_genetic)[1] <- "Genotype_1"
names(Lsmean_genetic)[2] <- "Genotype"
names(Lsmean_genetic)[4] <- "Species"

write.csv(Lsmean_genetic, file="Chapter2UP/Eighttrait_backlsmean_genetic group.csv",row.names = F)

##keep unique 
Lsmean_genetic_unique <- unique(Lsmean_genetic)
str(Lsmean_genetic_unique)
names(Lsmean_genetic_unique)



write.csv(Lsmean_genetic_unique, file="Chapter2UP/Eighttrait_backlsmean_genetic_unique.csv",row.names = F)

###find the first ten largeest genetic individuals in each traits 

str(Lsmean_genetic_unique)

Yld <- Lsmean_genetic_unique[,c(4,2:3,13,21)]

Yld_order <- Yld[order(Yld$Yld,decreasing=TRUE),]

write.csv(Yld_order,file="Chapter2UP/Yld_order.csv",row.names = F)

Bcirc <- Lsmean_genetic_unique[,c(4,2:3,6,14)]

Bcirc_order <- Bcirc[order(Bcirc$Bcirc,decreasing=TRUE),]

str(Bcirc_order)

write.csv(Bcirc_order,file="Chapter2UP/Bcirc_order.csv",row.names = F)



CCirc <- Lsmean_genetic_unique[,c(4,2:3,7,15)]

CCirc_order <- CCirc[order(CCirc$CCirc,decreasing=TRUE),]

str(CCirc_order)

write.csv(CCirc_order,file="Chapter2UP/CCirc_order.csv",row.names = F)

CmD_BI <- Lsmean_genetic_unique[,c(4,2:3,8,16)]

CmD_BI_order <- CmD_BI[order(CmD_BI$CmD_BI,decreasing=TRUE),]

str(CmD_BI_order)

write.csv(CmD_BI_order,file="Chapter2UP/CmD_BI_order.csv",row.names = F)

CmD_LI <- Lsmean_genetic_unique[,c(4,2:3,9,17)]

CmD_LI_order <- CmD_LI[order(CmD_LI$CmD_LI,decreasing=TRUE),]

str(CmD_LI_order)

write.csv(CmD_LI_order,file="Chapter2UP/CmD_LI_order.csv",row.names = F)


CmDW  <- Lsmean_genetic_unique[,c(4,2:3,10,18)]

CmDW_order <- CmDW[order(CmDW$CmDW,decreasing=TRUE),]

str(CmDW_order)

write.csv(CmDW_order,file="Chapter2UP/CmDW_order.csv",row.names = F)

Cml  <- Lsmean_genetic_unique[,c(4,2:3,11,19)]

Cml_order <- Cml[order(Cml$Cml,decreasing=TRUE),]

str(Cml_order)

write.csv(Cml_order,file="Chapter2UP/Cml_order.csv",row.names = F)

CmN  <- Lsmean_genetic_unique[,c(4,2:3,12,20)]

CmN_order <- CmN[order(CmN$CmN,decreasing=TRUE),]

str(CmN_order)

write.csv(CmN_order,file="Chapter2UP/CmN_order.csv",row.names = F)


#####get the lsmean for the different genotype  (169) for different traits (eight traits) 
##import the transformated datasets using boxcox method 

qualdat_check<- read.csv("data/traits1718normalited_Yld.csv")
str(qualdat_check)
qualdat.chap2ls <- qualdat_check[, -c(1:2,14:15,17:18)]
str(qualdat.chap2ls)
names(qualdat.chap2ls)

###import the dataset  with orignal souce 
Access.Source.traits <- read.csv("data/Accession.Source.traits_updated_edit.csv")
names(Access.Source.traits)

Access.Source.traitsls <- Access.Source.traits[,-c(1,3,17:30)]
names(Access.Source.traitsls)

###merge these two datasets using "Entry"

Acess.traitsls <- merge(Access.Source.traitsls,qualdat.chap2ls,by=c("Entry"),all=TRUE)
names(Acess.traitsls)

##rename the column names 

colnames(Acess.traitsls)[1] <- "Genotype"
colnames(Acess.traitsls)[5] <- "Genetic group"


###this one is the LS mean of different genotypic group for Figure S1
###this one is the LS mean of different genotypic group for Figure S1

###import the lsmean from SAS and Lsmean different for Figure S1
###import the lsmean from SAS and Lsmean different for Figure S1

#LsmeanSAS <- read.csv("Chapter2UP/Lsmean_up.csv")
LsmeanSAS <- read.csv("Chapter2UP/Lsmean_up-4.csv") #Lsmean_up-4 this one get from SAS 

head(LsmeanSAS) 
str(LsmeanSAS)
levels(as.factor(LsmeanSAS$filename))
LsmeanSAS$filename <- gsub('Lsmean.', '',LsmeanSAS$filename)
LsmeanSAS$filename <- gsub('.xlsx', '',LsmeanSAS$filename)
head(LsmeanSAS)
##back tansformation 

###import the lda

lda <- read.csv("data/lda.makeup_kg.csv")
lda_1 <-lda[-c(7,9,10,12,13),] 
str(lda_1)
##change the first column in to row.name

#rownames(lda_1) <- lda_1[,1]
#lda_1[,1] <- NULL
#lda_1

##add lda_1 to the datasets of smrylsmean

str(LsmeanSAS)

LsmeanSAS

lsmeanSASlda <- merge(LsmeanSAS, lda_1, by.x=c("filename"),by.y=c("X"), alll.x=TRUE)

##faction of the bcBack

bcBack <- function(x1, lda){
  ifelse(lda==0, exp(x1), (lda * x1 + 1)^(1/lda) )
}

str(lsmeanSASlda)

lsmeanSASlda$BackLSmeans <- bcBack(lsmeanSASlda$Estimate,lsmeanSASlda$lambda)

lsmeanSASlda$Backupper.CL <- bcBack(lsmeanSASlda$Upper, lsmeanSASlda$lambda)

lsmeanSASlda$Backlower.CL <- bcBack(lsmeanSASlda$Lower, lsmeanSASlda$lambda)


write.csv(lsmeanSASlda, file="Chapter2UP/lsmeanSASlda.csv",row.names = F)

####

str(lsmeanSASlda)

lsmeanSASlda[13:15] <- round(lsmeanSASlda[13:15],2)

str(lsmeanSASlda)

lsmeanSASlda$Lsmean_CI <- paste0(lsmeanSASlda$BackLSmeans," (",lsmeanSASlda$Backlower.CL,"-",lsmeanSASlda$Backupper.CL,")")

levels(as.factor(lsmeanSASlda$filename))

library(ggplot2)
levels(as.factor(lsmeanSASlda$DAPC_group))

for (k in levels(as.factor(lsmeanSASlda$filename))){
  data_1 = lsmeanSASlda[lsmeanSASlda$filename==k,]
  #data_1$variables <- factor(data_1$variables,levels=unique(data_1$variables[order(data_1$BackLSmeans)]))
  
  data_1$group_1 <- ifelse(data_1$DAPC_group == "M. xgiganteus 3x", "red", "blue")
  
  data_1$color <- ifelse(data_1$DAPC_group == "M. xgiganteus 3x", 0, 1)
  a <- ifelse(data_1$color == 0, "red", "blue")
  
  #data_1$group <- factor(data_1$group, levels = (Inds))
  
  plot <- ggplot(data = data_1,
                 aes(x= DAPC_group,y=BackLSmeans,label=Lsmean_CI)) +
    geom_point(aes(x= DAPC_group, y=BackLSmeans), stat="identity",colour="Blue", fill="white", alpha=3) +
    
    coord_cartesian(ylim=c(min(data_1$Backlower.CL),max(data_1$Backupper.CL)))+
    
    theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    
    geom_errorbar(aes(x=DAPC_group, ymin=Backlower.CL, ymax=Backupper.CL), width=0.2, colour=data_1$group_1, alpha=0.6, size=0.5)+
    geom_text(aes(label= Lsmean_CI),size=3,hjust=0.25,vjust=-0.5)+
    coord_flip(ylim=NULL)+
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major.x = element_line(colour = "lightgray", size = 0.5),
          panel.grid.minor.y = element_line(color = "lightgrey"),
          panel.border = element_rect(color = "lightgrey", fill = NA),
          axis.ticks.length.x = unit(0.50, "cm"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          axis.title.y=element_blank(),
          axis.text.y =element_text(size = 11,face = "bold"),
          axis.title=element_text(),
          axis.title.x=element_blank(),
          plot.caption=element_text(size=10,face="bold"),
          strip.text = element_text(face="bold", size=11),
          legend.title = element_text(size=10, face="bold"),
          legend.text = element_text(size=10, face="bold"),
          legend.position="bottom",axis.text.x=element_text(face="bold",size=9.0, angle=0))+
    labs(title= paste("Back-transformed Least-squared means and standard error for", k)) 
  ggsave(plot, path = "Chapter2UP/Mean_image_SAS_up", filename =paste(k, "pdf", sep = '.'),
         width=10,height = 6)
}


