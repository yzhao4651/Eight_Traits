###import all the .csv file from different file 
###emerging all the data together 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPNEW/")
#setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GPNEW/")
extension <- "csv"
fileNames <- Sys.glob(paste("GP*/Rsquare*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i], row.names=1)
  require(pastecs)
  options(scipen=100)
  options(digits=4)
  stat.desc(sample)
  mz1=stat.desc(sample)
  names <- rownames(mz1)
  rownames(mz1) <- NULL
  data <- cbind(names,mz1)
  mzList[[i]] = data.frame(data, filename = rep(fileNames[i], length(mz1)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
GP.Summ<- plyr::ldply(mzList, data.frame)
str(GP.Summ)
levels(GP.Summ$filename)
GP.Summ$filename <- gsub('.csv', '', GP.Summ$filename)
GP.Summ$filename <- gsub('Rsquare ', '', GP.Summ$filename)
GP.Summ$filename <- gsub('GP/', 'GPIM/', GP.Summ$filename)
GP.Summ$filename <- gsub('GP', '', GP.Summ$filename)
str(GP.Summ)
levels(as.factor(GP.Summ$filename))
##separated by the"/"
library(tidyverse)
GP.Summ.all <- separate(GP.Summ, filename,c("Ind.SNP","Trait"),"/")
GP.Summ.all$Trait <- gsub('GR', '', GP.Summ.all$Trait)
str(GP.Summ.all)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
write.csv(GP.Summ.all, file= "data/GPtotal3.Summ.all.csv")
levels(as.factor(GP.Summ.all$Trait))
levels(as.factor(GP.Summ.all$Ind.SNP))


GP.all.1 <- GP.Summ.all
levels(as.factor(GP.all.1$Trait))
levels(as.factor(GP.all.1$Ind.SNP))
levels(as.factor(GP.all.1$names))
GP.all.1 <- droplevels(GP.all.1)
GP.all.1$Ind.SNP <- as.factor(GP.all.1$Ind.SNP)

levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RR0"] <- "F116im/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo106R0"] <- "F106/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo116R0"] <- "F116/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo96R0"] <- "F96/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA102R0"] <- "OWA102/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA112R0"] <- "OWA112/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA122R0"] <- "OWA122/NOPC"
#levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA135R"] <- "OWA135/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC106R0"] <- "C106/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC116R0"] <- "C116/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC124R0"] <- "C124/NOPC"

levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RR1"] <- "F116im/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo106R1"] <- "F106/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo116R1"] <- "F116/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo96R1"] <- "F96/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA102R1"] <- "OWA102/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA112R1"] <- "OWA112/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA122R1"] <- "OWA122/PC1"
#levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA135R"] <- "OWA135/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC106R1"] <- "C106/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC116R1"] <- "C116/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC124R1"] <- "C124/PC1"

levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RR2"] <- "F116im/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo106R2"] <- "F106/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo116R2"] <- "F116/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo96R2"] <- "F96/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA102R2"] <- "OWA102/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA112R2"] <- "OWA112/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA122R2"] <- "OWA122/PC1~2"
#levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA135R"] <- "OWA135/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC106R2"] <- "C106/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC116R2"] <- "C116/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC124R2"] <- "C124/PC1~2"

levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RR3"] <- "F116im/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo106R3"] <- "F106/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo116R3"] <- "F116/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo96R3"] <- "F96/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA102R3"] <- "OWA102/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA112R3"] <- "OWA112/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA122R3"] <- "OWA122/PC1~3"
#levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA135R"] <- "OWA135/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC106R3"] <- "C106/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC116R3"] <- "C116/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC124R3"] <- "C124/PC1~3"


GP.all.2 <- separate(GP.all.1, Ind.SNP,c("Ind.SNP","PC"),"/")

levels(as.factor(GP.all.2$Trait))
levels(as.factor(GP.all.2$Ind.SNP))
levels(as.factor(GP.all.2$names))
levels(as.factor(GP.all.2$PC))
GP.all.2 <- droplevels(GP.all.2)

GP.all.2$Trait <- as.factor(GP.all.2$Trait)
write.csv(GP.all.2,file="data/GP.all.2.csv")
levels(as.factor(GP.all.2$Trait))

##get flowering traits
GP.all.flo <- GP.all.2[GP.all.2$Trait=="FD_1 "|GP.all.2$Trait=="FD_50. "|GP.all.2$Trait=="FM_1 "|GP.all.2$Trait=="FM_50 "|
                       GP.all.2$Trait== "fprind " |GP.all.2$Trait=="fprinGW "|GP.all.2$Trait=="fprinM " |GP.all.2$Trait=="fprinW "
                     |GP.all.2$Trait=="FW_1 "|GP.all.2$Trait=="FW_50 "|GP.all.2$Trait=="GFW_1 "
                     |GP.all.2$Trait=="GFW_50 "|GP.all.2$Trait=="GHW_1 "|GP.all.2$Trait=="GHW_50 "|GP.all.2$Trait=="HD_1 " |
                       GP.all.2$Trait=="HD_50. " 
                     |GP.all.2$Trait=="HM_1 " |GP.all.2$Trait=="HM_50 "|GP.all.2$Trait=="HW_1 " |GP.all.2$Trait=="HW_50 ",]
GP.all.flo <- droplevels(GP.all.flo)
str(GP.all.flo)
levels(as.factor(GP.all.flo$Trait))
GP.all.flo$Trait <- as.factor(GP.all.flo$Trait)
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="fprind "] <- "PC1/PC1.D"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="fprinGW "] <- "PC1/PC1.GW"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="fprinM "] <- "PC1/PC1.M"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="fprinW "] <- "PC1/PC1.W"

levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)== "FD_1 "] <- "FT1/FD_1"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="FM_1 "] <- "FT1/FM_1"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="GFW_1 "] <- "FT1/GFW_1"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="FW_1 "] <- "FT1/FW_1"

levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="FD_50. "] <- "FT50/FD_50"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="FM_50 "] <- "FT50/FM_50"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="FW_50 "] <- "FT50/FW_50"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="GFW_50 "] <- "FT50/GFW_50"

levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="HD_1 "] <- "HT1/HD_1"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)== "GHW_1 "] <- "HT1/GHW_1"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="HM_1 "] <- "HT1/HM_1"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="HW_1 "] <- "HT1/HW_1"

levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="GHW_50 "] <- "HT50/GHW_50"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="HD_50. "] <- "HT50/HD_50"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="HM_50 "] <- "HT50/HM_50"
levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="HW_50 "] <- "HT50/HW_50"
#levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="FD_50. "] <- "FD_50 "
#levels(GP.all.flo$Trait)[levels(GP.all.flo$Trait)=="HD_50. "] <- "HD_50 "
str(GP.all.flo)

GP.all.flo.2 <- separate(GP.all.flo, Trait,c("Trait","HFT"),"/")

GP.all.flo.2<- droplevels(GP.all.flo.2)
levels(as.factor(GP.all.flo.2$Trait))
levels(as.factor(GP.all.flo.2$Ind.SNP))
levels(as.factor(GP.all.flo.2$names))
levels(as.factor(GP.all.flo.2$HFT))
GP.all.flo.2 <- droplevels(GP.all.flo.2)
##write out the result
write.csv(GP.all.flo.2, file="data/GP.all.floup.csv")

### plot the image with FT_1
GP.all.flo.3 <- GP.all.flo.2[GP.all.flo.2$names=="mean",]

GP.all.flo.3 <- droplevels(GP.all.flo.3)
GP.all.flo.3$Ind.SNP<- factor(GP.all.flo.3$Ind.SNP,levels = c("F116im","F96","F106", "F116"))
GP.all.flo.2$HFT <- as.factor(GP.all.flo.2$HFT)
levels(GP.all.flo.2$HFT)
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="PC1.D"] <- "D"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="PC1.GW"] <- "GW"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="PC1.M"] <- "M"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="PC1.W"] <- "W"

levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)== "FD_1"] <- "D"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="FM_1"] <- "M"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="GFW_1"] <- "GW"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="FW_1"] <- "W"

levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="FD_50"] <- "D"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="FM_50"] <- "M"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="FW_50"] <- "W"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="GFW_50"] <- "GW"

levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="HD_1"] <- "D"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)== "GHW_1"] <- "GW"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="HM_1"] <- "M"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="HW_1"] <- "W"

levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="GHW_50"] <- "GW"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="HD_50"] <- "D"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="HM_50"] <- "M"
levels(GP.all.flo.2$HFT)[levels(GP.all.flo.2$HFT)=="HW_50"] <- "W"

levels(as.factor(GP.all.flo.2$HFT))
### plot the image with FT_1
GP.all.flo.4 <- GP.all.flo.2[GP.all.flo.2$names=="mean",]
GP.all.flo.4 <- droplevels(GP.all.flo.4)
str(GP.all.flo.4)
write.csv(GP.all.flo.4,file="Flotable/GP.mean.csv")
GP.all.flo.4$Ind.SNP<- factor(GP.all.flo.4$Ind.SNP,levels = c("F116im","F96", "F106", "F116"))
ALL <- ggplot(data=GP.all.flo.4, aes(x=PC, y=x, group=HFT, colour=HFT)) +
  geom_line() +
  geom_point()+facet_grid(Trait~Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
      axis.text=element_text(size = 12),
      axis.title=element_text(size = 12),
      axis.title.y=element_text(size=12),
      plot.caption=element_text(size=12),
      strip.text = element_text(face="bold", size=12),
      legend.title = element_text(size=12, face="bold"),
      legend.text = element_text(size=12, face="bold"),
      legend.position="bottom",
      axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies",title="The mean of prediction accuracies across SNPs combinations and Traits")

ALL.PC <- ggplot(data=GP.all.flo.4, aes(x=Ind.SNP, y=x, group=HFT, colour=HFT)) +
  geom_line() +
  geom_point()+facet_grid(Trait~PC)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies",title="The mean of prediction accuracies across SNPs combinations and Traits")

GP.all.flo.4 <- GP.all.flo.2[GP.all.flo.2$names=="std.dev",]
GP.all.flo.4 <- droplevels(GP.all.flo.4)
str(GP.all.flo.4)
write.csv(GP.all.flo.4,file="Flotable/GP.st.csv")
GP.all.flo.4$Ind.SNP<- factor(GP.all.flo.4$Ind.SNP,levels = c("F116im","F96", "F106", "F116"))
ALL.st <- ggplot(data=GP.all.flo.4, aes(x=PC, y=x, group=HFT, colour=HFT)) +
  geom_line() +
  geom_point()+facet_grid(Trait~Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Standard deviation of prediction accuracies",title="Standard deviation of prediction accuracies across SNPs combinations and Traits")


ALL.st.pc <- ggplot(data=GP.all.flo.4, aes(x=Ind.SNP, y=x, group=HFT, colour=HFT)) +
  geom_line() +
  geom_point()+facet_grid(Trait~PC)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Standard deviation of prediction accuracies",title="Standard deviation of prediction accuracies across SNPs combinations and Traits")

### plot the image with FT_1
GP.all.flo <- GP.all.flo[GP.all.flo$names=="mean",]
GP.all.flo.1 <- GP.all.flo[GP.all.flo$Trait=="FD_1 "| GP.all.flo$Trait=="FW_1 "| GP.all.flo$Trait=="GFW_1 "| GP.all.flo$Trait=="FM_1 ",]
GP.all.flo.1 <- droplevels(GP.all.flo.1)
str(GP.all.flo.1)
GP.all.flo.1$Ind.SNP<- factor(GP.all.flo.1$Ind.SNP,levels = c("F116im","F96", "F106", "F116"))
FT1 <- ggplot(data=GP.all.flo.1, aes(x=PC, y=x, group=Trait, colour=Trait)) +
  geom_line() +
  geom_point()+facet_grid(~Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies",title="Half Flowering time")

### plot the image with FT_50
GP.all.flo <- GP.all.flo[GP.all.flo$names=="mean",]
GP.all.flo.1 <- GP.all.flo[GP.all.flo$Trait=="FD_50 "| GP.all.flo$Trait=="FW_50 "| GP.all.flo$Trait=="GFW_50 "| GP.all.flo$Trait=="FM_50 ",]
GP.all.flo.1 <- droplevels(GP.all.flo.1)
str(GP.all.flo.1)
GP.all.flo.1$Ind.SNP<- factor(GP.all.flo.1$Ind.SNP,levels = c("F116im","F96", "F106", "F116"))
FT50 <- ggplot(data=GP.all.flo.1, aes(x=PC, y=x, group=Trait, colour=Trait)) +
  geom_line() +
  geom_point()+facet_grid(~Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies",title="Half Flowering time")


### plot the image with HT_1
GP.all.flo <- GP.all.flo[GP.all.flo$names=="mean",]
GP.all.flo.1 <- GP.all.flo[GP.all.flo$Trait=="HD_1 "| GP.all.flo$Trait=="HW_1 "| GP.all.flo$Trait=="GHW_1 "| GP.all.flo$Trait=="HM_1 ",]
GP.all.flo.1 <- droplevels(GP.all.flo.1)
str(GP.all.flo.1)
GP.all.flo.1$Ind.SNP<- factor(GP.all.flo.1$Ind.SNP,levels = c("F116im","F96", "F106", "F116"))
HT1 <- ggplot(data=GP.all.flo.1, aes(x=PC, y=x, group=Trait, colour=Trait)) +
  geom_line() +
  geom_point()+facet_grid(~Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies",title="Heading time")

### plot the image with HT_50
GP.all.flo <- GP.all.flo[GP.all.flo$names=="mean",]
GP.all.flo.1 <- GP.all.flo[GP.all.flo$Trait=="HD_50 "| GP.all.flo$Trait=="HW_50 "| GP.all.flo$Trait=="GHW_50 "| GP.all.flo$Trait=="HM_50 ",]
GP.all.flo.1 <- droplevels(GP.all.flo.1)
str(GP.all.flo.1)
GP.all.flo.1$Ind.SNP<- factor(GP.all.flo.1$Ind.SNP,levels = c("F116im","F96", "F106", "F116"))
HT50 <- ggplot(data=GP.all.flo.1, aes(x=PC, y=x, group=Trait, colour=Trait)) +
  geom_line() +
  geom_point()+facet_grid(~Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies",title="Half Heading time")

### plot the image with HT_50
GP.all.flo <- GP.all.flo[GP.all.flo$names=="mean",]
GP.all.flo.1 <- GP.all.flo[GP.all.flo$Trait=="PC1.D"| GP.all.flo$Trait=="PC1.W"| GP.all.flo$Trait=="PC1.GW"| GP.all.flo$Trait=="PC1.M",]
GP.all.flo.1 <- droplevels(GP.all.flo.1)
str(GP.all.flo.1)
GP.all.flo.1$Ind.SNP<- factor(GP.all.flo.1$Ind.SNP,levels = c("F116im","F96", "F106", "F116"))
PC <- ggplot(data=GP.all.flo.1, aes(x=PC, y=x, group=Trait, colour=Trait)) +
  geom_line() +
  geom_point()+facet_grid(~Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies",title="PC1 of Flowering time")

###combining all of the images together 

library(gridExtra)
library(ggpubr)
figall.5 <- ggarrange(HT1, HT50,FT1,FT50, PC, labels = c("A", "B","C","D","E","F"), ncol = 1, nrow = 6)


#t.13 <- GP.all.flo[GP.all.flo$Ind.SNP=="flo106"|GP.all.flo$Ind.SNP=="flo106R1"|
                     #GP.all.flo$Ind.SNP=="flo106R2"|GP.all.flo$Ind.SNP=="flo106R3"|GP.all.flo$Ind.SNP=="flo106R4"|
                     #GP.all.flo$Ind.SNP=="flo106R5",]
#t.13 <- t.13[t.13$Trait=="FD_1 ",]
#str(t.13)
#t.13 <- droplevels(t.13)
#str(t.13)
#t.13$Ind.SNP<- factor(t.13$Ind.SNP, levels= c("flo106", "flo106R1","flo106R2","flo106R3","flo106R4","flo106R5"))
#plot(t.13$Ind.SNP, t.13$x, type="line",col="blue", main="Mean of prediction accuracies of F106",ylab="Mean", xlab="Flo106",ylim=)

###get the OWA traits
write.csv(GP.all.2,file="data/GP.all.2.csv",row.names = 1)

GP.all.2 <- read.csv("data/GP.all.2.csv",row.names = 1)
str(GP.all.2)
##get flowering traits
levels(GP.all.2$Trait)
GP.all.OWA <- GP.all.2[GP.all.2$Trait=="OWA18 "|GP.all.2$Trait=="OWA19 ",]
GP.all.OWA <- droplevels(GP.all.OWA)
str(GP.all.OWA)
levels(as.factor(GP.all.OWA$Ind.SNP))
GP.all.OWA$Ind.SNP <- as.factor(GP.all.OWA$Ind.SNP)
levels(GP.all.OWA$Ind.SNP)[levels(GP.all.OWA$Ind.SNP)=="F116im"] <- "122im"
levels(GP.all.OWA$Ind.SNP)[levels(GP.all.OWA$Ind.SNP)=="OWA122"] <- "O122"
levels(GP.all.OWA$Ind.SNP)[levels(GP.all.OWA$Ind.SNP)=="OWA112"] <- "O112"
levels(GP.all.OWA$Ind.SNP)[levels(GP.all.OWA$Ind.SNP)=="OWA102"] <- "O102"
levels(as.factor(GP.all.OWA$Ind.SNP))
### 
GP.all.OWA.m <- GP.all.OWA[GP.all.OWA$names=="mean",]
str(GP.all.OWA.m)
write.csv(GP.all.OWA.m,file="OWA2/GP.chapter2.gp.mean.csv")
GP.all.OWA.m <- droplevels(GP.all.OWA.m)
GP.all.OWA.m$Ind.SNP<- factor(GP.all.OWA.m$Ind.SNP,levels = c("122im", "O122", "O112", "O102"))
str(GP.all.OWA.m)
library(ggplot2)

library(directlabels) ### this the package for the geom_dl function 
###this one is for PC were separated and Trait were group together 
OWA.Trait.group <- ggplot(data = GP.all.OWA.m, aes(x=Ind.SNP, y=x, group=Trait,colour=Trait)) +
  geom_line() +
  geom_point() + facet_grid(~PC) +
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = Trait,colour=Trait), method = list(dl.combine(method= "first.bumpup", method= "last.bumpup"), cex = 0.8))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies", title="The mean of prediction accuracies")
OWA.Trait.group

OWA.SNP.group <- ggplot(data = GP.all.OWA.m, aes(x=PC, y=x, group=Trait,colour=Trait)) +
  geom_line() +
  geom_point() + facet_grid(~Ind.SNP) +
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = Trait,colour=Trait), method = list(dl.combine(method= "first.bumpup", method= "last.bumpup"), cex = 0.8))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies", title="The mean of prediction accuracies")
OWA.SNP.group

###this one for stadard derivation 

GP.all.OWA.m <- GP.all.OWA[GP.all.OWA$names=="std.dev",]
str(GP.all.OWA.m)
write.csv(GP.all.OWA.m,file="OWA2/GP.chapter2.gp.std.csv")
GP.all.OWA.m <- droplevels(GP.all.OWA.m)
GP.all.OWA.m$Ind.SNP<- factor(GP.all.OWA.m$Ind.SNP,levels = c("122im", "O122", "O112", "O102"))
str(GP.all.OWA.m)
library(ggplot2)

library(directlabels) ### this the package for the geom_dl function 
###this one is for PC were separated and Trait were group together 
OWA.Trait.group <- ggplot(data = GP.all.OWA.m, aes(x=Ind.SNP, y=x, group=Trait,colour=Trait)) +
  geom_line() +
  geom_point() + facet_grid(~PC) +
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = Trait,colour=Trait), method = list(dl.combine(method= "first.bumpup", method= "last.bumpup"), cex = 0.8))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="The Standard deviation of prediction accuracies", title="The Standard deviation of prediction accuracies")
OWA.Trait.group

OWA.SNP.group <- ggplot(data = GP.all.OWA.m, aes(x=PC, y=x, group=Trait,colour=Trait)) +
  geom_line() +
  geom_point() + facet_grid(~Ind.SNP) +
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = Trait,colour=Trait), method = list(dl.combine(method= "first.bumpup", method= "last.bumpup"), cex = 0.8))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="The Standard deviation of prediction accuracies", title="Standard deviation of prediction accuracies")
OWA.SNP.group



###this one only for OWA with imputed SNPs
levels(GP.all.2$Trait)
GP.all.OWA <- GP.all.2[GP.all.2$Trait=="OWA ",]
GP.all.OWA <- droplevels(GP.all.OWA)
str(GP.all.OWA)
levels(as.factor(GP.all.OWA$Ind.SNP))
GP.all.OWA$Ind.SNP <- as.factor(GP.all.OWA$Ind.SNP)
levels(GP.all.OWA$Ind.SNP)[levels(GP.all.OWA$Ind.SNP)=="F116im"] <- "O122im"
levels(as.factor(GP.all.OWA$Ind.SNP))
### plot the image with FT_1
GP.all.OWA.m <- GP.all.OWA[GP.all.OWA$names=="mean",]
str(GP.all.OWA.m)
#write.csv(GP.all.OWA.m,file="OWA2/GP.chapter2.gp.mean.csv")
GP.all.OWA.m <- droplevels(GP.all.OWA.m)
#GP.all.OWA.m$Ind.SNP<- factor(GP.all.OWA.m$Ind.SNP,levels = c("O122im", "OWA122", "OWA112", "OWA102"))
str(GP.all.OWA.m)
library(ggplot2)

library(directlabels) ### this the package for the geom_dl function 
###this one is for PC were separated and Trait were group together 


OWA.Trait.group <- ggplot(data = GP.all.OWA.m, aes(x=PC, y=x,group=Trait)) +
  geom_line()+
  geom_point()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="Mean of prediction accuracies", title="The mean of prediction accuracies of Trait of OWA")
OWA.Trait.group

