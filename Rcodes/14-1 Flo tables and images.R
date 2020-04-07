##import the data file with total SNPs detected 
###making the table7 including all of the SNP combination with the hit flowering 
###read the dataset 
floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)
hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)
###get the table 
##combine the SNP and the Ind.SNP


alltable1 <- as.data.frame(table(floall.1[,c(6)]))
colnames(alltable1)[colnames(alltable1)=="Var1"] <- "Ind.SNP Combinations"
colnames(alltable1)[colnames(alltable1)=="Freq"] <- "SNP Totally detected"
hitflotable1 <- as.data.frame(table(hitflo[,c(6)]))
colnames(hitflotable1)[colnames(hitflotable1)=="Var1"] <- "Ind.SNP Combinations"
colnames(hitflotable1)[colnames(hitflotable1)=="Freq"] <- "SNPs Hit With Flower Genes"
##merge the two table 
t.2 <- plyr::join_all(list(alltable1, hitflotable1), by="Ind.SNP Combinations")
###wite out the table 
write.csv(t.2,file="Flotable/IndSNP.csv")


###imput the table and then making the bar images 

floall.1 <- read.csv("Flotraits2/floalltraits1.csv",row.names = 1)
floall.1 <- floall.1[!(floall.1$Ind.SNP=="F116+36088im"),]
floall.1 <- droplevels(floall.1)
str(floall.1)
#table.1 <- as.data.frame(table(floall.1$Ind.SNP))

floall.2 <-  unite_(floall.1, "SNP.PS", c("SNP","Ind.SNP"),sep=";")
floall.2.table <- as.data.frame(table(floall.2$SNP.PS))

t.2 <- floall.2.table[!(duplicated(floall.2.table$Var1)),]
###
t.1 <- separate(t.2, "Var1", c("SNP","Ind.SNP"),sep=";")

alltable1 <- as.data.frame(table(t.1[,c(2)]))
colnames(alltable1)[colnames(alltable1)=="Var1"] <- "Non.missing"
colnames(alltable1)[colnames(alltable1)=="Freq"] <- "SNPs Totally Detected"

##get the hit with flowering 

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv",row.names = 1)
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hit.2 <-  unite_(hitflo, "SNP.PS", c("SNP","Ind.SNP"),sep=";")
hit.2.table <- as.data.frame(table(hit.2$SNP.PS))

t.3 <- hit.2.table[!(duplicated(hit.2.table$Var1)),]
###
t.4 <- separate(t.3, "Var1", c("SNP","Ind.SNP"),sep=";")

hitflotable1 <- as.data.frame(table(t.4[,c(2)]))
colnames(hitflotable1)[colnames(hitflotable1)=="Var1"] <- "Non.missing"
colnames(hitflotable1)[colnames(hitflotable1)=="Freq"] <- "SNPs Hit With Flower Genes"
##merge the two table 
t.5 <- plyr::join_all(list(alltable1, hitflotable1), by="Non.missing")
str(t.5)
library(reshape2)
library(reshape)
t.16<- melt(t.5, id.vars = c("Non.missing"))
str(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With Flower Genes", "SNPs Totally Detected"))
write.csv(t.16, file="data/SNPs.M.csv")
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
t.16[t.16 == 0] <- NA
t.16$Non.missing <- factor(t.16$Non.missing,levels = c("F96+4814", "F106+4185", "F116+3098"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
NON.mis <- ggplot(data=t.16, aes(x=Non.missing, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+
  geom_text(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) +
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
  labs(y="The number of QTNs",title="The non-missing combinations")

###
######this table made by hnad "Ind.SNP.Combinations"
SNPcombs <- read.csv("Flotable/Ind.SNP.Combinations1.csv")
str(SNPcombs)
t.15 <- SNPcombs[4:6]
colnames(t.15)[colnames(t.15)=="SNP.Totally.Detected.1"] <- "SNP.Totally.Detected"
colnames(t.15)[colnames(t.15)=="SNPs.Hit.With.Flower.Genes.1"] <- "SNPs.Hit.With.Flower.Genes"
library(reshape2)
library(reshape)
t.18 <- melt(t.15, id.vars = c("Ind.SNP.Combinations"))
str(t.18)
library(ggplot2)
t.18$variable<- factor(t.18$variable, levels= c("SNPs.Hit.With.Flower.Genes", "SNP.Totally.Detected"))

#t.18$Method <- factor(t.18$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
t.18[t.18 == 0] <- NA
t.18$Ind.SNP.Combinations <- factor(t.18$Ind.SNP.Combinations,levels = c("Imputed", "Non-missing", "Imputed and non-missing"))
#t.18$Trait.name<- factor(t.18$Trait.name,levels = c("D", "W", "GW", "M"))

IM.NON.Over<- ggplot(data=t.18, aes(x=Ind.SNP.Combinations, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+
  geom_text(aes(label=t.18$value),size=4, position = position_stack(vjust = 0.5)) +
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
  labs(y="The number of QTNs",title="The Imputed and Non-missing Combinations")

###write this images into one images 

library(ggpubr)
#theme_set(theme_pubr())
#install.packages("gridExtra")
library(gridExtra)

fig.all <- ggarrange(IM.NON.Over, NON.mis,labels = c("A", "B"), ncol = 2, nrow =1 )

##working on the PS_N and PS_Y, and overlap
##working on the PS_N and PS_Y, and overlap
##working on the PS_N and PS_Y, and overlap
##working on the PS_N and PS_Y, and overlap
##import the data set
floall.1 <- read.csv("Flotraits/floalltraits1.csv")
str(floall.1)
levels
##select the method from mrMLM 
floall.1 <- floall.1[floall.1$Software=="mrMLMM",]
floall.1 <- droplevels(floall.1)
str(floall.1)
hitflo.1 <- read.csv("Flotraits/hitflowithotherspeciesup.csv")

###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

###this one for getting the PS_Y,PS_N, and PS_Y and PS_N 
floall.2 <-  unite_(floall.1, "SNP.PS", c("SNP","PS"),sep=";")
floall.2.table <- as.data.frame(table(floall.2$SNP.PS))
floall.2.table.1 <- floall.2.table[1]
###
t.1 <- separate(floall.2.table.1, "Var1", c("SNP","PS"),sep=";")
##get all of the duplicated 
t.2 <- t.1[duplicated(t.1$SNP),]

write.csv(t.2,file="Flotable/t.2.csv")
t.3 <- t.1[t.1$SNP %in% t.2$SNP,]
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
write.csv(t.3,file="Flotable/t.3.csv")
###change the PS_N and PS_Y
t.3$PS <- as.factor(t.3$PS)
levels(t.3$PS)[levels(t.3$PS)=="PS_Y"] <- "PS_Y and PS_N"
levels(t.3$PS)[levels(t.3$PS)=="PS_N"] <- "PS_Y and PS_N"
###remove the duplicated value 
t.4<- unique(t.3)
###join it together with the 
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
t.flo <- rbind(t.4,t.5)

####this one working for the hit with flo
##get the PS
hitflo.2 <-  unite_(hitflo, "SNP.PS", c("SNP","PS"),sep=";")
hitflo.2.table <- as.data.frame(table(hitflo.2$SNP.PS))
hitflo.2.table.1 <- hitflo.2.table[1]
###
t.6 <- separate(hitflo.2.table.1, "Var1", c("SNP","PS"),sep=";")
##get all of the duplicated 
t.7 <- t.6[duplicated(t.6$SNP),]
t.8 <- t.6[t.6$SNP %in% t.7$SNP,]
###change the PS_N and PS_Y
t.8$PS <- as.factor(t.8$PS)
levels(t.8$PS)[levels(t.8$PS)=="PS_Y"] <- "PS_Y and PS_N"
levels(t.8$PS)[levels(t.8$PS)=="PS_N"] <- "PS_Y and PS_N"
###remove the duplicated value 
t.10 <- unique(t.8)
###join it together with the 
t.9 <- t.6[!(t.6$SNP %in% t.7$SNP),]
t.hit <- rbind(t.9,t.10)

###get the table 
alltable1 <- as.data.frame(table(t.flo$PS))
colnames(alltable1)[colnames(alltable1)=="Var1"] <- "PS"
colnames(alltable1)[colnames(alltable1)=="Freq"] <- "SNP Totally Detected"
hitflotable1 <- as.data.frame(table(t.hit$PS))
colnames(hitflotable1)[colnames(hitflotable1)=="Var1"] <- "PS"
colnames(hitflotable1)[colnames(hitflotable1)=="Freq"] <- "SNPs Hit With Flower Genes"
##merge the two table 
t.2 <- plyr::join_all(list(alltable1, hitflotable1), by="PS")
###wite out the table 
write.csv(t.2,file="Flotable/PS.csv")
###imput the table and then making the bar images 
###this table made by hnad "Ind.SNP.Combinations"
str(t.2)
library(reshape2)
library(reshape)
t.16<- melt(t.2, id.vars = c("PS"))
str(t.16)
#levels(t.16$PS)[levels(t.16$PS)=="PS_N;PS_Y"] <- "PS_N and PS_Y"
str(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With Flower Genes", "SNP Totally Detected"))
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
t.16[t.16 == 0] <- NA
t.16$PS <- factor(t.16$PS,levels = c("PS_N", "PS_Y", "PS_Y and PS_N"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))

PS <- ggplot(data=t.16, aes(x=t.16$PS, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="PS_N:No Population Structure\nPS_Y:Population Structure\nPS_Y and PS_N:Overlap",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12), axis.text.y=element_text(size=12))+
  labs(x="PS_N=Population Structure not in Model,PS_Y=Population Structure in Model,PS_Y and PS_N=Overlap",y="The number of QTNs",title="Population structure with or without")

####Method and SNP.combinations 
####Method and SNP.combinations 
####Method and SNP.combinations 
filename <- read.csv("Flotraits/floalltraits2.csv")
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="F116im+36088"] <- "F116+36088im"
levels(filename$Method)
str(filename)
str(filename)
###change all of traits to the four temporal scales 
levels(filename$Trait.name)[levels(filename$Trait.name)=="HD_1"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HW_1"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="GHW_1"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HM_1"] <- "M"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HD_50"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HW_50"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="GHW_50"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HM_50"] <- "M"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FD_1"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FW_1"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="GFW_1"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FM_1"] <- "M"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FD_50"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FW_50"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="GFW_50"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FM_50"] <- "M"
levels(filename$Trait.name)[levels(filename$Trait.name)=="fprind"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="fprinW"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="fprinGW"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="fprinM"] <- "M"
library(tidyr)
library(dplyr)
##select the column need for study
t <- as.data.frame(ftable(filename[,c(7,10,11)]))
t.T.1 <- droplevels(t)
str(t.T.1)

###unit the SNP. method and trait.nmes 
t.T.1.u <-  unite_(t.T.1, "IndMT", c("Ind.SNP","Method", "Trait.name"),sep=";")
##change the column name Freq to  "SNPs.Total.detected"
colnames(t.T.1.u)[colnames(t.T.1.u)=="Freq"] <- "SNPs.Total.detected"
###import the SNP hit with flowering genes 
hitflowering <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
##select the SNPs hit with flowering but in the totole detected SNPs. 
hitflo <- filename[filename$SNP %in% hitflowering$SNP,]
hitflo <- droplevels(hitflo)
tflo <- as.data.frame(ftable(hitflo[,c(7,10,11)]))
nlevels(hitflo$SNP)
###unite the Method and Trait.name
###unit the SNP. method and trait.nmes in order to use to merge with the total detected SNPs
t.T.2.u <-  unite_(tflo, "IndMT", c("Ind.SNP", "Method", "Trait.name"),sep=";")
## change the column name Freq to "SNPs.hit.with.flowering.genes"
colnames(t.T.2.u)[colnames(t.T.2.u)=="Freq"] <- "SNPs.hit.with.flowering.genes"
##merge two files t.T.1.u, and t.T.2.u
#t.1 <- merge(x=t.T.1.u, y=t.T.2.u, by ="IndMT", by.x=TRUE)## this one does not work 
t.2 <- plyr::join_all(list(t.T.1.u,t.T.2.u), by="IndMT")
t.2[is.na(t.2)] <- 0
##separate MT to method and Trait.name after merge
t.13 <- separate(t.2, "IndMT", c("Ind.SNP","Method", "Trait.name"),sep=";")
###write out the table 
write.csv(t.13, file="Flotable/Trait.Med.SNP.csv")
##transform from wide form to vertical form. 
library(reshape2)
library(reshape)
t.14 <- melt(t.13, id.vars = c("Ind.SNP","Method","Trait.name"))
str(t.14)
###F116+36088im
levels(as.factor(t.14$Method))
t.15 <- t.14[t.14$Method=="MLM+SUPER" | t.14$Method=="CMLM+SUPER",]
t.16 <- t.15[t.15$Ind.SNP=="F116+36088im",]
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.detected"))
t.16$Ind.SNP <- factor(t.16$Ind.SNP,levels = c("F116+36088im", "F106+4185", "F116+3077", "F122+2272"))
t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
total.1 <- ggplot(data=t.16, aes(x=Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) + facet_grid(~Method+Ind.SNP)+
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
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="F116+36088im")

###106+4185, F116+3077, and F122+2772
t.14 <- melt(t.13, id.vars = c("Ind.SNP","Method","Trait.name"))
levels(as.factor(t.14$Method))
t.15 <- t.14[t.14$Method=="MLM+SUPER" | t.14$Method=="CMLM+SUPER",]
t.16 <- t.15[!(t.15$Ind.SNP=="F116+36088im"),]
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.detected"))
t.16$Ind.SNP <- factor(t.16$Ind.SNP,levels = c("F96+4814","F106+4185", "F116+3098"))
t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
MLMSUPER.CMLMSUPER<- ggplot(data=t.16, aes(x=Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) + facet_grid(~Method+Ind.SNP)+
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
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="MLM+SUPER, CMLM+SUPER")

###all Multiple methods and other combinations. 
###for all of multiple methods
t.14 <- melt(t.13, id.vars = c("Ind.SNP","Method","Trait.name"))
str(t.14)
t.16 <- t.14[!(t.14$Method=="CMLM+SUPER"|t.14$Method=="MLM+SUPER"),]
#levels(as.factor(t.15$Method))
#t.16 <- t.15[t.15$Ind.SNP=="F122+2272",]
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.detected"))
write.csv(t.16, file="data/SNPs.M.csv")
t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
t.16[t.16 == 0] <- NA
t.16$Ind.SNP <- factor(t.16$Ind.SNP,levels = c("F116+36088im", "F96+4814", "F106+4185","F116+3098" ))
t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
allmultiple <- ggplot(data=t.16, aes(x=Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) + facet_grid(Method~Ind.SNP)+
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
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="all")

###getting the multiple methods 

floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)
floallmethod <-  unite_(floall.1, "ThreeMeth", c("SUPER","FarmCPU", "mrMLM.1","mrMLM.2"),sep=" ")
floallmethod.table <- as.data.frame(table(floallmethod$ThreeMeth))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("SUPER","FarmCPU", "mrMLM.1","mrMLM.2"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflomethod$ThreeMeth))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.2 <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.2[is.na(t.2)] <- 0
write.csv(t.2, file="Flotable/allmethod.csv")
##separate MT to method and Trait.name after merge
#t.13 <- separate(t.2, "IndMT", c("Ind.SNP","Method", "Trait.name"),sep=";")


###getting the method FarmCPU with SUPER
floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)
floallmethod <-  unite_(floall.1, "ThreeMeth", c("SUPER","FarmCPU"),sep=" ")
floallmethod.table <- as.data.frame(table(floallmethod$ThreeMeth))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("SUPER","FarmCPU"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflomethod$ThreeMeth))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.2.1 <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.2.1[is.na(t.2.1)] <- 0
write.csv(t.2.1, file="Flotable/SUPERFarmCPU.csv")

####SUPER with mrMLM.1
floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)
floallmethod <-  unite_(floall.1, "ThreeMeth", c("SUPER","mrMLM.1"),sep=" ")
floallmethod.table <- as.data.frame(table(floallmethod$ThreeMeth))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("SUPER","mrMLM.1"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflomethod$ThreeMeth))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.2.2 <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.2.2[is.na(t.2.2)] <- 0
write.csv(t.2.2, file="Flotable/SUPERmrMLM.1.csv")

####SUPER with mrMLM.2
floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)
floallmethod <-  unite_(floall.1, "ThreeMeth", c("SUPER","mrMLM.2"),sep=" ")
floallmethod.table <- as.data.frame(table(floallmethod$ThreeMeth))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("SUPER","mrMLM.2"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflomethod$ThreeMeth))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.2.3 <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.2.3[is.na(t.2.3)] <- 0
write.csv(t.2.3, file="Flotable/SUPERmrMLM.2.csv")

####FarmCPU with mrMLM.1
floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)
floallmethod <-  unite_(floall.1, "ThreeMeth", c("FarmCPU","mrMLM.1"),sep=" ")
floallmethod.table <- as.data.frame(table(floallmethod$ThreeMeth))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("FarmCPU","mrMLM.1"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflomethod$ThreeMeth))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.2.4 <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.2.4[is.na(t.2.4)] <- 0
write.csv(t.2.4, file="Flotable/FarmCPUmrMLM.1.csv")


####FarmCPU with mrMLM.2
floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)
floallmethod <-  unite_(floall.1, "ThreeMeth", c("FarmCPU","mrMLM.2"),sep=" ")
floallmethod.table <- as.data.frame(table(floallmethod$ThreeMeth))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("FarmCPU","mrMLM.2"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflomethod$ThreeMeth))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.2.5 <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.2.5[is.na(t.2.5)] <- 0
write.csv(t.2.5, file="Flotable/FarmCPUmrMLM.2.csv")

####mrMLM.1, mrMLM.2

floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)
floallmethod <-  unite_(floall.1, "ThreeMeth", c("mrMLM.1","mrMLM.2"),sep=" ")
floallmethod.table <- as.data.frame(table(floallmethod$ThreeMeth))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("mrMLM.1","mrMLM.2"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflomethod$ThreeMeth))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.2.6 <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.2.6[is.na(t.2.6)] <- 0
write.csv(t.2.6, file="Flotable/mrMLMmrMLM.2.csv")

###combina all of them together 

T.method <- rbind(t.2.1,t.2.2,t.2.3,t.2.4,t.2.5,t.2.6)

write.csv(T.method, file="Flotable/T.method.csv")



####SUPER
####SUPER
floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)
floallmethod.table <- as.data.frame(table(floall.1$SUPER))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"
hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)
hitflomethod.table <- as.data.frame(table(hitflo$SUPER))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"
t.super <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.super[is.na(t.super)] <- 0
write.csv(t.super, file="Flotable/SUPER.csv")

####FarmCPU

floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)

floallmethod.table <- as.data.frame(table(floall.1$FarmCPU))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("mrMLM.1","mrMLM.2"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflo$FarmCPU))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.Farm <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.Farm[is.na(t.Farm)] <- 0
write.csv(t.Farm, file="Flotable/FarmCPU.csv")

####mrMLM.1

floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)

floallmethod.table <- as.data.frame(table(floall.1$mrMLM.1))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("mrMLM.1","mrMLM.1"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflo$mrMLM.1))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.mrMLM.1 <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.mrMLM.1[is.na(t.mrMLM.1)] <- 0
write.csv(t.mrMLM.1, file="Flotable/mrMLM.1.csv")

####mrMLM.2
floall.1 <- read.csv("Flotraits2/flounique3table.csv")
str(floall.1)

floallmethod.table <- as.data.frame(table(floall.1$mrMLM.2))
colnames(floallmethod.table)[colnames(floallmethod.table)=="Freq"] <- "SNPs.Total.Detected"

hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("mrMLM.1","mrMLM.2"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflo$mrMLM.2))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.hit.with.flowering.genes"

t.mrMLM.2 <- plyr::join_all(list(floallmethod.table,hitflomethod.table), by="Var1")
t.mrMLM.2[is.na(t.mrMLM.2)] <- 0
write.csv(t.mrMLM.2, file="Flotable/mrMLM.2.csv")

###combina all of them together 

T.eachmethod <- rbind(t.super,t.Farm,t.mrMLM.1,t.mrMLM.2)

write.csv(T.eachmethod, file="Flotable/each.method.csv")

###making the table for each pair of the methods 
##import the data 
pairmethod <- read.csv("Flotable/T.methodUP.csv")
str(pairmethod)
###for all of multiple methods
t.14 <- melt(pairmethod, id.vars = c("Var1"))
str(t.14)
t.16 <- t.14
#levels(as.factor(t.15$Method))
#t.16 <- t.15[t.15$Ind.SNP=="F122+2272",]
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.Detected"))
#write.csv(t.16, file="data/SNPs.M.csv")
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
t.16[t.16 == 0] <- NA
t.16$vaue <- factor(t.16$value,order(t.16$value))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
Paired.Meth<- ggplot(data=t.16, aes(x=reorder(Var1,-value), y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) +
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
  labs(y="The number of QTNs",title="The number of SNPs from any two methods")

###avoid the labels overlapped change geom_text_repel to  geom_text_repel
t.16 <- t.14
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.Detected"))
#write.csv(t.16, file="data/SNPs.M.csv")
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))

t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.Detected"))
t.16$vaue <- factor(t.16$value,order(t.16$value))
library(ggrepel)
Paired.Meth.1 <- ggplot(data=t.16, aes(x=reorder(Var1,-value), y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text_repel(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) +
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
  labs(y="The number of QTNs",title="The number of SNPs from any two methods")

###calculate the PC1 and with their temperal scales 
###PC.1 and Day
###PC1.D and Day
###
##import the data set
floall.1 <- read.csv("Flotraits/floalltraits1.csv")
str(floall.1)
levels(floall.1$Trait.name)
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="HD_1"] <- "D"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="HW_1"] <- "W"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="GHW_1"] <- "GW"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="HM_1"] <- "M"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="HD_50"] <- "D"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="HW_50"] <- "W"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="GHW_50"] <- "GW"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="HM_50"] <- "M"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="FD_1"] <- "D"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="FW_1"] <- "W"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="GFW_1"] <- "GW"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="FM_1"] <- "M"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="FD_50"] <- "D"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="FW_50"] <- "W"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="GFW_50"] <- "GW"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="FM_50"] <- "M"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="fprind"] <- "PC1.D"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="fprinW"] <- "PC1.W"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="fprinGW"] <- "PC1.GW"
levels(floall.1$Trait.name)[levels(floall.1$Trait.name)=="fprinM"] <- "PC1.M"
levels(floall.1$Trait.name)
##select the Day and PC1.D 
floall.2 <- floall.1[floall.1$Trait.name=="D"|floall.1$Trait.name=="PC1.D",]
floall.2 <- droplevels(floall.2)
str(floall.2)

floall.2 <-  unite_(floall.2, "SNP.PS", c("SNP","Trait.name"),sep=";")
floall.2.table <- as.data.frame(table(floall.2$SNP.PS))
floall.2.table.1 <- floall.2.table[1]
###
t.1 <- separate(floall.2.table.1, "Var1", c("SNP","Trait.name"),sep=";")
##get all of the duplicated 
t.2 <- t.1[duplicated(t.1$SNP),]

#write.csv(t.2,file="Flotable/t.2.csv")
t.3 <- t.1[t.1$SNP %in% t.2$SNP,]
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
#write.csv(t.3,file="Flotable/t.3.csv")
###change the D and PC1_D
t.3$Trait.name <- as.factor(t.3$Trait.name)
levels(t.3$Trait.name)[levels(t.3$Trait.name)=="PC1.D"] <- "PC1.D and D"
levels(t.3$Trait.name)[levels(t.3$Trait.name)=="D"] <- "PC1.D and D"
###remove the duplicated value 
t.4<- unique(t.3)
###join it together with the 
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
t.flo.D <- rbind(t.4,t.5)

####this one working for the hit with flo
##get the PS
hitflo.1 <- read.csv("Flotraits/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
t.hit.D <- t.flo.D[t.flo.D$SNP %in% hitflo.1$SNP,]
str(t.hit.D)


#### Week and PC1.week 
#### Week and PC1.week
#### Week and PC1.week
floall.2 <- floall.1[floall.1$Trait.name=="W"|floall.1$Trait.name=="PC1.W",]
floall.2 <- droplevels(floall.2)
str(floall.2)

floall.2 <-  unite_(floall.2, "SNP.PS", c("SNP","Trait.name"),sep=";")
floall.2.table <- as.data.frame(table(floall.2$SNP.PS))
floall.2.table.1 <- floall.2.table[1]
###
t.1 <- separate(floall.2.table.1, "Var1", c("SNP","Trait.name"),sep=";")
##get all of the duplicated 
t.2 <- t.1[duplicated(t.1$SNP),]

#write.csv(t.2,file="Flotable/t.2.csv")
t.3 <- t.1[t.1$SNP %in% t.2$SNP,]
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
#write.csv(t.3,file="Flotable/t.3.csv")
###change the W and PC1_W
t.3$Trait.name <- as.factor(t.3$Trait.name)
levels(t.3$Trait.name)[levels(t.3$Trait.name)=="PC1.W"] <- "PC1.W and W"
levels(t.3$Trait.name)[levels(t.3$Trait.name)=="W"] <- "PC1.W and W"
###remove the duplicated value 
t.4<- unique(t.3)
###join it together with the 
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
t.flo.W <- rbind(t.4,t.5)

####this one working for the hit with flo
##get the PS
hitflo.1 <- read.csv("Flotraits/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
t.hit.W  <- t.flo.W[t.flo.W$SNP %in% hitflo.1$SNP,]



#### GWeek and PC1.week 
#### GWeek and PC1.week 
#### GWeek and PC1.week 
floall.2 <- floall.1[floall.1$Trait.name=="GW"|floall.1$Trait.name=="PC1.GW",]
floall.2 <- droplevels(floall.2)
str(floall.2)

floall.2 <-  unite_(floall.2, "SNP.PS", c("SNP","Trait.name"),sep=";")
floall.2.table <- as.data.frame(table(floall.2$SNP.PS))
floall.2.table.1 <- floall.2.table[1]
###
t.1 <- separate(floall.2.table.1, "Var1", c("SNP","Trait.name"),sep=";")
##get all of the duplicated 
t.2 <- t.1[duplicated(t.1$SNP),]

#write.csv(t.2,file="Flotable/t.2.csv")
t.3 <- t.1[t.1$SNP %in% t.2$SNP,]
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
#write.csv(t.3,file="Flotable/t.3.csv")
###change the GW and PC1_GW
t.3$Trait.name <- as.factor(t.3$Trait.name)
levels(t.3$Trait.name)[levels(t.3$Trait.name)=="PC1.GW"] <- "PC1.GW and GW"
levels(t.3$Trait.name)[levels(t.3$Trait.name)=="GW"] <- "PC1.GW and GW"
###remove the duplicated value 
t.4<- unique(t.3)
###join it together with the 
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
t.flo.GW <- rbind(t.4,t.5)

####this one working for the hit with flo
##get the PS
hitflo.1 <- read.csv("Flotraits/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
t.hit.GW <- t.flo.GW[t.flo.GW$SNP %in% hitflo.1$SNP,]
str(t.hit.GW)


#### Month and PC1.M 
#### Month and PC1.M 
#### Month and PC1.M 
floall.2 <- floall.1[floall.1$Trait.name=="M"|floall.1$Trait.name=="PC1.M",]
floall.2 <- droplevels(floall.2)
str(floall.2)
levels(floall.2$Trait.name)

floall.2 <-  unite_(floall.2, "SNP.PS", c("SNP","Trait.name"),sep=";")
floall.2.table <- as.data.frame(table(floall.2$SNP.PS))
floall.2.table.1 <- floall.2.table[1]
###
t.1 <- separate(floall.2.table.1, "Var1", c("SNP","Trait.name"),sep=";")
##get all of the duplicated 
t.2 <- t.1[duplicated(t.1$SNP),]

#write.csv(t.2,file="Flotable/t.2.csv")
t.3 <- t.1[t.1$SNP %in% t.2$SNP,]
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
#write.csv(t.3,file="Flotable/t.3.csv")
###change the M and PC1_M
t.3$Trait.name <- as.factor(t.3$Trait.name)
levels(t.3$Trait.name)[levels(t.3$Trait.name)=="PC1.M"] <- "PC1.M and M"
levels(t.3$Trait.name)[levels(t.3$Trait.name)=="M"] <- "PC1.M and M"
###remove the duplicated value 
t.4<- unique(t.3)
###join it together with the 
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
t.flo.M <- rbind(t.4,t.5)

####this one working for the hit with flo
##get the PS
hitflo.1 <- read.csv("Flotraits/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
t.hit.M <- t.flo.M[t.flo.M$SNP %in% hitflo.1$SNP,]
str(t.hit.M)


t.hit <- rbind(t.hit.D,t.hit.W,t.hit.GW,t.hit.M)
t.flo <- rbind(t.flo.D,t.flo.W,t.flo.GW,t.flo.M)


###get the table 
alltable1 <- as.data.frame(table(t.flo$Trait.name))
colnames(alltable1)[colnames(alltable1)=="Var1"] <- "Trait.name"
colnames(alltable1)[colnames(alltable1)=="Freq"] <- "SNP Totally Detected"
hitflotable1 <- as.data.frame(table(t.hit$Trait.name))
colnames(hitflotable1)[colnames(hitflotable1)=="Var1"] <- "Trait.name"
colnames(hitflotable1)[colnames(hitflotable1)=="Freq"] <- "SNPs Hit With Flower Genes"
##merge the two table 
t.2 <- plyr::join_all(list(alltable1, hitflotable1), by="Trait.name")
###wite out the table 
write.csv(t.2,file="Flotable/Trait.name.csv")
###imput the table and then making the bar images 
###this table made by hnad "Ind.SNP.Combinations"
str(t.2)
library(reshape2)
library(reshape)
t.16<- melt(t.2, id.vars = c("Trait.name"))
str(t.16)
#levels(t.16$PS)[levels(t.16$PS)=="PS_N;PS_Y"] <- "PS_N and PS_Y"
str(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With Flower Genes", "SNP Totally Detected"))
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
#t.16[t.16 == 0] <- NA
t.16$Trait.name <- factor(t.16$Trait.name,levels = c("D", "PC1.D", "PC1.D and D",
                                     "W", "PC1.W", "PC1.W and W",
                                     "GW", "PC1.GW", "PC1.GW and GW",
                                     "M", "PC1.M", "PC1.M and M"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
###using \n to change another row for text 
PC <- ggplot(data=t.16, aes(x=t.16$Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month\nPC1:First Principal Component",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D:Day  W:Week  GW:Group Week  M:Month ",y="The number of QTNs",title="PC1 and four temporal scales")

###get the total number of the SNPs in each of the temporal scales 
###Heading_1
###Heading_1
###Heading_1
filename <- read.csv("Flotraits2/flounique3tableforalltypetemporalscale.csv")
str(filename)
hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- filename[filename$SNP %in% hitflo.1$SNP,]
str(hitflo)
names(filename)

alltable1.H1 <- as.data.frame(table(filename$Heading_1))
colnames(alltable1.H1)[colnames(alltable1.H1)=="Var1"] <- "Type"
colnames(alltable1.H1)[colnames(alltable1.H1)=="Freq"] <- "SNP Totally Detected"
alltable1.H1$Trait.name <- paste("Heading_1")
hitflotable1.H1 <- as.data.frame(table(hitflo$Heading_1))
colnames(hitflotable1.H1)[colnames(hitflotable1.H1)=="Var1"] <- "Type"
colnames(hitflotable1.H1)[colnames(hitflotable1.H1)=="Freq"] <- "SNPs Hit With Flower Genes"
hitflotable1.H1$Trait.name <- paste("Heading_1")
###Heading_50
###Heading_50
###Heading_50
filename <- read.csv("Flotraits2/flounique3tableforalltypetemporalscale.csv")
str(filename)
hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- filename[filename$SNP %in% hitflo.1$SNP,]
str(hitflo)
names(filename)

alltable1.H5 <- as.data.frame(table(filename$Heading_50))
colnames(alltable1.H5)[colnames(alltable1.H5)=="Var1"] <- "Type"
colnames(alltable1.H5)[colnames(alltable1.H5)=="Freq"] <- "SNP Totally Detected"
alltable1.H5$Trait.name <- paste("Heading_50")
hitflotable1.H5 <- as.data.frame(table(hitflo$Heading_50))
colnames(hitflotable1.H5)[colnames(hitflotable1.H5)=="Var1"] <- "Type"
colnames(hitflotable1.H5)[colnames(hitflotable1.H5)=="Freq"] <- "SNPs Hit With Flower Genes"
hitflotable1.H5$Trait.name <- paste("Heading_50")
###Flowering_50
###Flowering_1
###Flowering_1
filename <- read.csv("Flotraits2/flounique3tableforalltypetemporalscale.csv")
str(filename)
hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- filename[filename$SNP %in% hitflo.1$SNP,]
str(hitflo)
names(filename)

alltable1.F1 <- as.data.frame(table(filename$Flowering_1))
colnames(alltable1.F1)[colnames(alltable1.F1)=="Var1"] <- "Type"
colnames(alltable1.F1)[colnames(alltable1.F1)=="Freq"] <- "SNP Totally Detected"
alltable1.F1$Trait.name <- paste("Flowering_1")
hitflotable1.F1 <- as.data.frame(table(hitflo$Flowering_1))
colnames(hitflotable1.F1)[colnames(hitflotable1.F1)=="Var1"] <- "Type"
colnames(hitflotable1.F1)[colnames(hitflotable1.F1)=="Freq"] <- "SNPs Hit With Flower Genes"
hitflotable1.F1$Trait.name <- paste("Flowering_1")

###Flowering_50
###Flowering_50
###Flowering_50
filename <- read.csv("Flotraits2/flounique3tableforalltypetemporalscale.csv")
str(filename)
hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- filename[filename$SNP %in% hitflo.1$SNP,]
str(hitflo)
names(filename)

alltable1.F5 <- as.data.frame(table(filename$Flowering_50))
colnames(alltable1.F5)[colnames(alltable1.F5)=="Var1"] <- "Type"
colnames(alltable1.F5)[colnames(alltable1.F5)=="Freq"] <- "SNP Totally Detected"
alltable1.F5$Trait.name <- paste("Flowering_50")
hitflotable1.F5 <- as.data.frame(table(hitflo$Flowering_50))
colnames(hitflotable1.F5)[colnames(hitflotable1.F5)=="Var1"] <- "Type"
colnames(hitflotable1.F5)[colnames(hitflotable1.F5)=="Freq"] <- "SNPs Hit With Flower Genes"
hitflotable1.F5$Trait.name <- paste("Flowering_50")

###PCA
###PCA
###PCA
filename <- read.csv("Flotraits2/flounique3tableforalltypetemporalscale.csv")
str(filename)
hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- filename[filename$SNP %in% hitflo.1$SNP,]
str(hitflo)
names(filename)

alltable1.PCA <- as.data.frame(table(filename$PCA))
colnames(alltable1.PCA)[colnames(alltable1.PCA)=="Var1"] <- "Type"
colnames(alltable1.PCA)[colnames(alltable1.PCA)=="Freq"] <- "SNP Totally Detected"
alltable1.PCA$Trait.name <- paste("PCA")
hitflotable1.PCA <- as.data.frame(table(hitflo$PCA))
colnames(hitflotable1.PCA)[colnames(hitflotable1.PCA)=="Var1"] <- "Type"
colnames(hitflotable1.PCA)[colnames(hitflotable1.PCA)=="Freq"] <- "SNPs Hit With Flower Genes"
hitflotable1.PCA$Trait.name <- paste("PCA")

#####combina 

hit <- rbind(hitflotable1.H1,hitflotable1.H5,hitflotable1.F1,hitflotable1.F5,hitflotable1.PCA)
hit.unit <- unite_(hit,"all",c("Type","Trait.name"),sep="/")
flo <- rbind(alltable1.H1,alltable1.H5,alltable1.F1,alltable1.F5,alltable1.PCA)
flo.unit <- unite_(flo,"all",c("Type","Trait.name"),sep="/")
all <- plyr::join_all(list(flo.unit,hit.unit),by="all")
allsep <- separate(all,"all",c("Type","Trait.name"),sep="/")
allsep.1 <- allsep[!(allsep$Type==""),]

###making the shape
t.14 <- melt(allsep.1,id.vars=c("Type","Trait.name"))

t.16 <- t.14[t.14$Trait.name=="Heading_1",]
t.16 <- droplevels(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With Flower Genes", "SNP Totally Detected"))
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
#t.16[t.16 == 0] <- NA
#t.16$Trait.name <- factor(t.16$Trait.name,levels = c("Heading_1", "Heading_50",
                                                     #"Flowering_1", "Flowering_50","PCA"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
###using \n to change another row for text 
H1 <- ggplot(data=t.16, aes(x=t.16$Type, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text_repel(aes(label=t.16$value),size=4,position = position_stack(vjust = 0.5))+ 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D:Day  W:Week  GW:Group Week  M:Month ",y="The number of QTNs",title="Heading time")

####Half heading time 
t.14 <- melt(allsep.1,id.vars=c("Type","Trait.name"))

t.16 <- t.14[t.14$Trait.name=="Heading_50",]
t.16 <- droplevels(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With Flower Genes", "SNP Totally Detected"))
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
#t.16[t.16 == 0] <- NA
#t.16$Trait.name <- factor(t.16$Trait.name,levels = c("Heading_1", "Heading_50",
#"Flowering_1", "Flowering_50","PCA"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
###using \n to change another row for text 
H5 <- ggplot(data=t.16, aes(x=t.16$Type, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text_repel(aes(label=t.16$value),size=4,position = position_stack(vjust = 0.5))+ 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D:Day  W:Week  GW:Group Week  M:Month ",y="The number of QTNs",title="Half Heading time")

###First flowering time 

t.16 <- t.14[t.14$Trait.name=="Flowering_1",]
t.16 <- droplevels(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With Flower Genes", "SNP Totally Detected"))
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
#t.16[t.16 == 0] <- NA
#t.16$Trait.name <- factor(t.16$Trait.name,levels = c("Heading_1", "Heading_50",
#"Flowering_1", "Flowering_50","PCA"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
###using \n to change another row for text 
F1 <- ggplot(data=t.16, aes(x=t.16$Type, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text_repel(aes(label=t.16$value),size=4,position = position_stack(vjust = 0.5))+ 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D:Day  W:Week  GW:Group Week  M:Month ",y="The number of QTNs",title="First Flowering time")

###Half Flowering time  

t.16 <- t.14[t.14$Trait.name=="Flowering_50",]
t.16 <- droplevels(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With Flower Genes", "SNP Totally Detected"))
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
#t.16[t.16 == 0] <- NA
#t.16$Trait.name <- factor(t.16$Trait.name,levels = c("Heading_1", "Heading_50",
#"Flowering_1", "Flowering_50","PCA"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
###using \n to change another row for text 
F5 <- ggplot(data=t.16, aes(x=t.16$Type, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text_repel(aes(label=t.16$value),size=4,position = position_stack(vjust = 0.5))+ 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D:Day  W:Week  GW:Group Week  M:Month ",y="The number of QTNs",title="Half Flowering time")

###PCA
###PCA

t.16 <- t.14[t.14$Trait.name=="PCA",]
t.16 <- droplevels(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With Flower Genes", "SNP Totally Detected"))
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
#t.16[t.16 == 0] <- NA
#t.16$Trait.name <- factor(t.16$Trait.name,levels = c("Heading_1", "Heading_50",
#"Flowering_1", "Flowering_50","PCA"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
###using \n to change another row for text 
PC1 <- ggplot(data=t.16, aes(x=t.16$Type, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text_repel(aes(label=t.16$value),size=4,position = position_stack(vjust = 0.5))+ 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D:Day  W:Week  GW:Group Week  M:Month ",y="The number of QTNs",title="First Principle Component")


####calculating all of all types 


###Heading_1
###Heading_1
###Heading_1
filename <- read.csv("Flotraits2/flounique3tableforalltypetemporalscale.csv")
str(filename)
filename <- filename[,c(1,7:11)]
filename.1 <- melt(filename,id.vars=c("SNP"))
###unite all together 
levels(filename.1$Trait.name)
levels(filename.1$value)[levels(filename.1$value)=="W;D"] <- "D;W"
levels(filename.1$value)[levels(filename.1$value)=="W;W"] <- "W"
levels(filename.1$value)[levels(filename.1$value)=="W;D;M"] <- "D;W;M" 
levels(filename.1$value)[levels(filename.1$value)=="D;W;W"] <- "D;W"
levels(filename.1$value)[levels(filename.1$value)=="W;W;M"] <- "W;M"
levels(filename.1$value)[levels(filename.1$value)=="D;W;W;M"] <- "D;W;M"
levels(filename.1$value)[levels(filename.1$value)=="D;W;GW"] <- "D;GW;W"
levels(filename.1$value)[levels(filename.1$value)=="W;D;GW"] <- "D;GW;W"
levels(filename.1$value)[levels(filename.1$value)=="D;W;GW;M"] <- "D;GW;W;M"
levels(filename.1$value)[levels(filename.1$value)== "W;D;GW;M"] <- "D;GW;W;M"
levels(filename.1$value)[levels(filename.1$value)=="W;GW;M"] <- "GW;W;M"
levels(filename.1$value)[levels(filename.1$value)=="W;GW"] <- "GW;W"


hitflo.1 <- read.csv("Flotraits2/hitflowithotherspeciesup.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- filename[filename$SNP %in% hitflo.1$SNP,]
str(hitflo)
hitflo.1 <- melt(hitflo,id.vars=c("SNP"))


levels(hitflo.1$value)
levels(hitflo.1$value)[levels(hitflo.1$value)=="W;D"] <- "D;W"
levels(hitflo.1$value)[levels(hitflo.1$value)=="W;W"] <- "W"
levels(hitflo.1$value)[levels(hitflo.1$value)=="W;D;M"] <- "D;W;M" 
levels(hitflo.1$value)[levels(hitflo.1$value)=="D;W;W"] <- "D;W"
levels(hitflo.1$value)[levels(hitflo.1$value)=="W;W;M"] <- "W;M"
levels(hitflo.1$value)[levels(hitflo.1$value)=="D;W;W;M"] <- "D;W;M"
levels(hitflo.1$value)[levels(hitflo.1$value)=="D;W;GW"] <- "D;GW;W"
levels(hitflo.1$value)[levels(hitflo.1$value)=="W;D;GW"] <- "D;GW;W"
levels(hitflo.1$value)[levels(hitflo.1$value)=="D;W;GW;M"] <- "D;GW;W;M"
levels(hitflo.1$value)[levels(hitflo.1$value)== "W;D;GW;M"] <- "D;GW;W;M"
levels(hitflo.1$value)[levels(hitflo.1$value)=="W;GW;M"] <- "GW;W;M"
levels(hitflo.1$value)[levels(hitflo.1$value)=="W;GW"] <- "GW;W"


alltable1.all <- as.data.frame(table(filename.1$value))
colnames(alltable1.all)[colnames(alltable1.all)=="Var1"] <- "Type"
colnames(alltable1.all)[colnames(alltable1.all)=="Freq"] <- "SNP Totally Detected"

hitflotable1.all <- as.data.frame(table(hitflo.1$value))
colnames(hitflotable1.all)[colnames(hitflotable1.all)=="Var1"] <- "Type"
colnames(hitflotable1.all)[colnames(hitflotable1.all)=="Freq"] <- "SNPs Hit With Flower Genes"

#####combine 


all.t <- plyr::join_all(list(alltable1.all,hitflotable1.all),by="Type")
write.csv(all.t, file="Flotable/temporalscalefromalltype.csv")
allsep.1 <- all.t[!(all.t$Type==""),]

all.t.1 <- melt(allsep.1,id.vars=c("Type"))

t.16 <- droplevels(all.t.1)
levels(t.16$Type)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With Flower Genes", "SNP Totally Detected"))
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
#t.16[t.16 == 0] <- NA
#t.16$Trait.name <- factor(t.16$Trait.name,levels = c("Heading_1", "Heading_50",
#"Flowering_1", "Flowering_50","PCA"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
###using \n to change another row for text 
T <- ggplot(data=t.16, aes(x=t.16$Type, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4,position = position_stack(vjust = 0.5))+ 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.position="bottom",
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D:Day  W:Week  GW:Group Week  M:Month ",y="The number of QTNs",title="Types from all")


###get total type of four type of temporal scale

###write all of the image into one 

library(ggpubr)
#theme_set(theme_pubr())
#install.packages("gridExtra")
library(gridExtra)
##this one is not very good.
#Temptype1 <- ggarrange(H1,H5,F1,F5,PC1,T,labels = c("A", "B","c","D","E","T"), ncol = 2, nrow =3, common.legend = TRUE, legend="bottom")



