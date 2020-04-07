
###making the image with different methods and traits and ind.SNP
###three ways with method, ind.SNP, and traits
filename <- read.csv("Charpter2/all.chapter2.csv")
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="143+36088im"] <- "124+36088im"
levels(filename$Trait.name)

levels(filename$Trait.name)[levels(filename$Trait.name)=="Bcirc_cm"] <- "Bcirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CCirc_cm"] <- "CCirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_BI_mm"] <- "CmD_BI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_LI_mm"] <- "CmD_LI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmDW_g"] <- "CmDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Cml_cm"] <- "Cml"
levels(filename$Trait.name)[levels(filename$Trait.name)=="SDW_kg"] <- "SDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Yld_kg"] <- "Yld"

str(filename)
t.H.1 <- as.data.frame(ftable(filename[,c(12,15,16)]))
t.H.1.nozero <- t.H.1[!(t.H.1$Freq==0),]
write.csv(t.H.1,file="Charpter2/table.Ind.SNP.Methodzero.csv")
write.csv(t.H.1.nozero,file="Charpter2/table.Ind.SNP.Method.csv")


###
t.H.1.MCM <- t.H.1[t.H.1$Method=="CMLM+SUPER" | t.H.1$Method=="MLM+SUPER" | t.H.1$Method=="GLM" | t.H.1$Method=="rrBLUP",]
#t.H.1.MCM <- t.H.1.MCM[t.H.1.MCM$Ind.SNP=="F116+36088im",]
t.H.1.MCM$Freq[t.H.1.MCM$Freq==0] <- NA
t.H.1.MCM <- droplevels(t.H.1.MCM)
t.H.1.MCM
str(t.H.1.MCM)
t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("124+36088im", "C106+4202", "C116+3293", "C125+2562"))
INdmethod <- ggplot(data=t.H.1.MCM, aes(x=Trait.name, y=Freq, fill=Trait.name)) +
  geom_bar(stat="identity",width=0.75, position=position_dodge(width=0.8),
           linetype = 0)+ facet_grid(Method~Ind.SNP)+
  geom_text(aes(x=t.H.1.MCM$Trait.name, y=t.H.1.MCM$Freq,label=t.H.1.MCM$Freq),size=4)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12), 
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=10,angle=45), axis.text.y=element_text(size=12))+ 
       labs(y="The number of QTNs")

t.H.1.MCM <- t.H.1[!(t.H.1$Method=="CMLM+SUPER" | t.H.1$Method=="MLM+SUPER" | t.H.1$Method=="GLM" | t.H.1$Method=="rrBLUP"),]
#t.H.1.MCM <- t.H.1.MCM[t.H.1.MCM$Ind.SNP=="F116+36088im",]
t.H.1.MCM$Freq[t.H.1.MCM$Freq==0] <- NA
t.H.1.MCM <- droplevels(t.H.1.MCM)
t.H.1.MCM
str(t.H.1.MCM)
t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("124+36088im", "C106+4202", "C116+3293", "C125+2562"))
INdMulmethod <- ggplot(data=t.H.1.MCM, aes(x=Trait.name, y=Freq, fill=Trait.name)) +
  geom_bar(stat="identity",width=0.75, position=position_dodge(width=0.8),
           linetype = 0)+ facet_grid(Method~Ind.SNP)+
  geom_text(aes(x=t.H.1.MCM$Trait.name, y=t.H.1.MCM$Freq,label=t.H.1.MCM$Freq),size=4)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=10,angle=45), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs")

###no method only ind.SNP

filename <- read.csv("data/all.chapter2.csv")
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="143+36088im"] <- "124+36088im"
levels(filename$Trait.name)
levels(filename$Trait.name)[levels(filename$Trait.name)=="Bcirc_cm"] <- "Bcirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CCirc_cm"] <- "CCirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_BI_mm"] <- "CmD_BI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_LI_mm"] <- "CmD_LI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmDW_g"] <- "CmDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Cml_cm"] <- "Cml"
levels(filename$Trait.name)[levels(filename$Trait.name)=="SDW_kg"] <- "SDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Yld_kg"] <- "Yld"
str(filename)
t.H.1.MCM <- as.data.frame(table(filename[,c(12,16)]))
##write out the result
write.csv(t.H.1.MCM,file="Charpter2/table.Ind.SNP.csv")


#t.H.1.MCM <- t.H.1[t.H.1$Method=="CMLM+SUPER" | t.H.1$Method=="MLM+SUPER" | t.H.1$Method=="GLM" | t.H.1$Method=="rrBLUP",]
#t.H.1.MCM <- t.H.1.MCM[t.H.1.MCM$Ind.SNP=="F116+36088im",]
t.H.1.MCM$Freq[t.H.1.MCM$Freq==0] <- NA
t.H.1.MCM <- droplevels(t.H.1.MCM)
t.H.1.MCM
str(t.H.1.MCM)
t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("124+36088im", "C106+4202", "C116+3293", "C125+2562"))
Trait.x <- ggplot(data=t.H.1.MCM, aes(x=Trait.name, y=Freq, fill=Trait.name)) +
  geom_bar(stat="identity",width=0.75, position=position_dodge(width=0.8),
           linetype = 0)+ facet_grid(~Ind.SNP)+
  geom_text(aes(x=t.H.1.MCM$Trait.name, y=t.H.1.MCM$Freq,label=t.H.1.MCM$Freq),size=4)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=10,angle=45), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs")

###trying to change the x-axis
#t.H.1.MCM <- t.H.1[t.H.1$Method=="CMLM+SUPER" | t.H.1$Method=="MLM+SUPER" | t.H.1$Method=="GLM" | t.H.1$Method=="rrBLUP",]
#t.H.1.MCM <- t.H.1.MCM[t.H.1.MCM$Ind.SNP=="F116+36088im",]
t.H.1.MCM$Freq[t.H.1.MCM$Freq==0] <- NA
t.H.1.MCM <- droplevels(t.H.1.MCM)
t.H.1.MCM
str(t.H.1.MCM)
t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("124+36088im", "C106+4202", "C116+3293", "C125+2562"))
Trait.x <- ggplot(data=t.H.1.MCM, aes(x=Trait.name, y=Freq, fill=Trait.name)) +
  geom_bar(stat="identity",width=0.75, position=position_dodge(width=0.8),
           linetype = 0)+ facet_grid(~Ind.SNP)+
  geom_text(aes(x=t.H.1.MCM$Trait.name, y=t.H.1.MCM$Freq,label=t.H.1.MCM$Freq),size=4)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=10,angle=45), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs")

##Remove the SNPs from GLM and rrBLUP for CmN trait because they did not include the K matrix 
## draw the line for the trend of the number of the SNP-traits associations instead of using the bar
install.packages("directlabels")
library(directlabels)
library(ggplot2)
library(directlabels)
filename <- read.csv("Charpter2/all.chapter2.csv")
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="143+36088im"] <- "124+36088im"
levels(filename$Trait.name)
#filename <- filename[!(filename$Method=="GLM"|filename=="rrBLUP"),]
#filename <- droplevels(filename)
t.H.1.MCM <- as.data.frame(table(filename[,c(12,16)]))
SNP.x <- ggplot(data=t.H.1.MCM, aes(x=Ind.SNP, y=Freq, group=Trait.name, colour=Trait.name)) +
  geom_line() +
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = Trait.name,colour=Trait.name), method = list(dl.combine(method= "first.bumpup", method= "last.bumpup"), cex = 0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=0), axis.text.y=element_text(size=12))+
  labs(y="The number of QTNs",title="The number of QTNs across the combinations of SNPs and accessions")

###the SNPs-traits associations from Method

filename <- read.csv("Charpter2/all.chapter2.csv")
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="143+36088im"] <- "124+36088im"
str(filename)
levels(filename$Trait.name)
#filename <- filename[!(filename$Method=="GLM"|filename=="rrBLUP"),]
#filename <- droplevels(filename)
t.H.1.M <- as.data.frame(table(filename[,c(15:16)]))
#t.H.1.M$Freq[t.H.1.M$Freq==0] <- NA
Methods.x <- ggplot(data=t.H.1.M, aes(x=Method, y=Freq, group=Trait.name, colour=Trait.name)) +
  geom_line() +
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = Trait.name,colour=Trait.name), method = list(dl.combine(method= "first.bumpup", method= "last.bumpup"), cex = 0.8)) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=90), axis.text.y=element_text(size=12))+
  labs(y="The number of QTNs",title="The number of QTNs from all methods")

##draw the bar diagram 
###the SNPs-traits associations from Method
filename <- read.csv("Charpter2/all.chapter2.csv")
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="143+36088im"] <- "124+36088im"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Bcirc_cm"] <- "Bcirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CCirc_cm"] <- "CCirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_BI_mm"] <- "CmD_BI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_LI_mm"] <- "CmD_LI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmDW_g"] <- "CmDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Cml_cm"] <- "Cml"
levels(filename$Trait.name)[levels(filename$Trait.name)=="SDW_kg"] <- "SDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Yld_kg"] <- "Yld"
str(filename)
levels(filename$Trait.name)
#filename <- filename[!(filename$Method=="GLM"|filename=="rrBLUP"),]
#filename <- droplevels(filename)
t.H.1.MCM <- as.data.frame(table(filename[,c(15:16)]))
##write out the dataset 
write.csv(t.H.1.MCM,file="Charpter2/table.method.csv")
t.H.1.MCM$Freq[t.H.1.MCM$Freq==0] <- NA
t.H.1.MCM <- droplevels(t.H.1.MCM)
t.H.1.MCM
str(t.H.1.MCM)
#t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("124+36088im", "C106+4202", "C116+3293", "C125+2562"))
Method.bar.x <- ggplot(data=t.H.1.MCM, aes(x=Trait.name, y=Freq, fill=Trait.name)) +
  geom_bar(stat="identity",width=0.75, position=position_dodge(width=0.8),
           linetype = 0)+ facet_grid(~Method)+
  geom_text(aes(x=t.H.1.MCM$Trait.name, y=t.H.1.MCM$Freq,label=t.H.1.MCM$Freq),size=4)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=10,angle=90), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs",title="The number of QTNs from all methods")



###no method only ind.SNP
##only have 36088im
filename <- read.csv("data/all.chapter2.csv")
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="143+36088im"] <- "124+36088im"
levels(filename$Trait.name)
levels(filename$Trait.name)[levels(filename$Trait.name)=="Bcirc_cm"] <- "Bcirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CCirc_cm"] <- "CCirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_BI_mm"] <- "CmD_BI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_LI_mm"] <- "CmD_LI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmDW_g"] <- "CmDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Cml_cm"] <- "Cml"
levels(filename$Trait.name)[levels(filename$Trait.name)=="SDW_kg"] <- "SDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Yld_kg"] <- "Yld"
str(filename)
t.H.1 <- as.data.frame(table(filename[,c(12,16)]))
t.H.1.MCM <- t.H.1[t.H.1$Ind.SNP=="124+36088im",]
#t.H.1.MCM <- t.H.1.MCM[t.H.1.MCM$Ind.SNP=="F116+36088im",]
t.H.1.MCM$Freq[t.H.1.MCM$Freq==0] <- NA
t.H.1.MCM <- droplevels(t.H.1.MCM)
t.H.1.MCM
str(t.H.1.MCM)
#t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("124+36088im", "C106+4202", "C116+3293", "C125+2562"))
Trait.x <- ggplot(data=t.H.1.MCM, aes(x=Trait.name, y=Freq, fill=Trait.name)) +
  geom_bar(stat="identity",width=0.75, position=position_dodge(width=0.8),
           linetype = 0)+ facet_grid(~Ind.SNP)+
  geom_text(aes(x=t.H.1.MCM$Trait.name, y=t.H.1.MCM$Freq,label=t.H.1.MCM$Freq),size=4)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=12,angle=0), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs")


Trait.x.1 <- ggplot(data=t.H.1.MCM, aes(x=Ind.SNP, y=Freq, fill=Trait.name)) +
  geom_bar(stat="identity",width=0.75, position=position_dodge(width=0.8),
           linetype = 0)+ facet_grid(~Trait.name)+
  geom_text(aes(x=t.H.1.MCM$Ind.SNP, y=t.H.1.MCM$Freq,label=t.H.1.MCM$Freq),size=4)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=12,angle=0), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs")

###only have non-missing SNPs
t.H.1 <- as.data.frame(table(filename[,c(12,16)]))
t.H.1.MCM <- t.H.1[!(t.H.1$Ind.SNP=="124+36088im"),]
#t.H.1.MCM <- t.H.1.MCM[t.H.1.MCM$Ind.SNP=="F116+36088im",]
t.H.1.MCM$Freq[t.H.1.MCM$Freq==0] <- NA
t.H.1.MCM <- droplevels(t.H.1.MCM)
t.H.1.MCM
str(t.H.1.MCM)
t.H.1.MCM$Freq[t.H.1.MCM$Freq==0] <- NA
t.H.1.MCM <- droplevels(t.H.1.MCM)
t.H.1.MCM
str(t.H.1.MCM)
t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("C106+4202", "C116+3293", "C125+2562"))
INdx <- ggplot(data=t.H.1.MCM, aes(x=Ind.SNP, y=Freq, fill=Ind.SNP)) +
  geom_bar(stat="identity",width=0.75, position=position_dodge(width=0.8),
           linetype = 0)+ facet_grid(~Trait.name)+
  geom_text(aes(x=t.H.1.MCM$Ind.SNP, y=t.H.1.MCM$Freq,label=t.H.1.MCM$Freq),size=4)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=10,angle=45), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs")


###change traits to be X.

t.H.1.MCM <- t.H.1[!(t.H.1$Method=="CMLM+SUPER" | t.H.1$Method=="MLM+SUPER" | t.H.1$Method=="GLM" | t.H.1$Method=="rrBLUP"),]
#t.H.1.MCM <- t.H.1.MCM[t.H.1.MCM$Ind.SNP=="F116+36088im",]
t.H.1.MCM$Freq[t.H.1.MCM$Freq==0] <- NA
t.H.1.MCM <- droplevels(t.H.1.MCM)
t.H.1.MCM
str(t.H.1.MCM)
t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("124+36088im", "C106+4202", "C116+3293", "C125+2562"))
INdMulmethod <- ggplot(data=t.H.1.MCM, aes(x=Trait.name, y=Freq, fill=Trait.name)) +
  geom_bar(stat="identity",width=0.75, position=position_dodge(width=0.8),
           linetype = 0)+ facet_grid(Method~Ind.SNP)+
  geom_text(aes(x=t.H.1.MCM$Trait.name, y=t.H.1.MCM$Freq,label=t.H.1.MCM$Freq),size=4)+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=10,angle=45), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs")



###calculate the total number of SNPs from different SNPs and overlaped by two SNPs
str(filename)
SNP.trait <- as.data.frame(table(filename[,c(1,16)]))
SNP.trait.nozero <- SNP.trait[!(SNP.trait$Freq==0),]
###write out the table 
write.csv(SNP.trait.nozero, file="Charpter2/SNP.trait.nozero.csv")### and then got the paired number through the pivot table in excel

##import the pair numbers of two traits 

paired <- read.csv("Charpter2/SNP.paired.traits.csv")
paired$group <- paired$Trait1
levels(paired$group)
str(paired)
paired.u <-  unite_(paired, "Two.Traits", c("Trait1", "Trait2"),sep="+")
paired.u.charp <-  paired.u[!(paired.u$group=="Flo"| paired.u$group=="OWA"| paired.u$group=="SRD"),]
str(paired.u.charp)
##reshape the dataset from wide to tall
paired.u.charp <- paired.u.charp[order(paired.u.charp$overlap),]
paired.u.charp$Two.Traits <- factor(paired.u.charp$Two.Traits,levels=paired.u.charp$Two.Traits[order(paired.u.charp$overlap)])
x$name <- factor(x$name, levels = x$name[order(x$val)])
paired.charper2 <- ggplot(data=paired.u.charp, aes(x=Two.Traits, y=overlap,fill=overlap)) +
  geom_bar(stat="identity",width=0.5,color="blue")+
  geom_text(aes(label=paired.u.charp$overlap),size=4, position = position_stack(vjust = 0.5), color="white") + 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=90), axis.text.y=element_text(size=12))+ labs(y="The number of Overlaped QTNs")

###import the correlation coeffient 

paired <- read.csv("Charpter2/Paired.coefficient.csv")
str(paired)
paired$group <- paired$Trait1
levels(paired$group)
str(paired)
paired.u <-  unite_(paired, "Two.Traits", c("Trait1", "Trait2"),sep="+")
#paired.u.charp <-  paired.u[!(paired.u$group=="Flo"| paired.u$group=="OWA"| paired.u$group=="SRD"),]
str(paired.u)
##reshape the dataset from wide to tall
paired.u <- paired.u[order(paired.u$COA),]
paired.u$Two.Traits <- factor(paired.u$Two.Traits,levels=paired.u$Two.Traits[order(paired.u$COA)])
#x$name <- factor(x$name, levels = x$name[order(x$val)])
paired.charper2.coa <- ggplot(data=paired.u, aes(x=Two.Traits, y=COA,fill=COA)) +
  geom_bar(stat="identity",width=0.5,color="blue")+
  geom_text(aes(label=paired.u$COA),size=4, position = position_stack(vjust = 0.5), color="white") + 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=90), axis.text.y=element_text(size=12))+ labs(y="The correlation of coefficient")

###image COA and overlapSNP together 
paired.COA <- read.csv("Charpter2/Paired.coefficient.csv")
paired.u.COA <-  unite_(paired.COA, "Two.Traits", c("Trait1", "Trait2"),sep="+")
paired.u.COA <- paired.u.COA[order(paired.u.COA$Two.Traits),]
paired.u.COA$Two.Traits
paired.SNP <- read.csv("Charpter2/SNP.paired.traits.1.csv")
paired.u.SNP <- paired.u.SNP[order(paired.u.SNP$Two.Traits),]
paired.u.SNP <-  unite_(paired.SNP, "Two.Traits", c("Trait1", "Trait2"),sep="+")
paired.u.SNP <- paired.u.SNP[order(paired.u.SNP$Two.Traits),]
paired.u.SNP.sort.names <- paired.u.SNP[order(paired.u.SNP$overlap),]
paired.u.SNP$Two.Traits
paired.u.SNP$Two.Traits %in% paired.u.COA$Two.Traits
##merger both 
Pair.COA.SNP.1 <- merge(paired.u.SNP, paired.u.COA, all = FALSE, all.y= TRUE, by.y="Two.Traits")
str(Pair.COA.SNP.1)

library(reshape2)
library(reshape)
Pair.COA.SNP <- melt(Pair.COA.SNP.1, id.vars = c("Two.Traits"))
str(Pair.COA.SNP)
levels(Pair.COA.SNP$variable)[levels(Pair.COA.SNP$variable)=="overlap"] <- "overlapped SNPs"
levels(Pair.COA.SNP$variable)[levels(Pair.COA.SNP$variable)=="COA"] <- "Correlation of coefficient"

Pair.COA.SNP$Two.Traits <- factor(Pair.COA.SNP$Two.Traits,levels=Pair.COA.SNP$Two.Traits[order(Pair.COA.SNP$COA)])
Pair.COA.SNP$variable<- factor(Pair.COA.SNP$variable, levels= c("Correlation of coefficient","overlapped SNPs"))
Pair.COA.SNP$Two.Traits <- factor(Pair.COA.SNP$Two.Traits,levels=Pair.COA.SNP$Two.Traits[order(paired.u.SNP$overlap)])
write.csv(Pair.COA.SNP,file="Charpter2/Pair.COA.SNP.csv")

COA.SNP <- ggplot(data=Pair.COA.SNP, aes(x=Two.Traits, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(labels=c("Correlation of coefficient","overlapped SNPs"))+   
  geom_text(aes(label=Pair.COA.SNP$value),size=4, position = position_stack(vjust = 0.5)) + 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=90), axis.text.y=element_text(size=12))+
  labs(y="The Correlation coefficient & the number of Overlapped SNP",title="The Correlation coefficient & the number of Overlapped SNP of Nine traits")






