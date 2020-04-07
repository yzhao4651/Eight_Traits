##import the data file with total SNPs detected 

###making the table7 including all of the SNP combination with the hit OWA 
###read the dataset 
OWAall.1 <- read.csv("OWA2/OWAuniquetable.csv")
str(OWAall.1)
hitOWA.1 <- read.csv("OWA2/AHG.owa.hit.T.csv")
###find the SNP hit with OWA genes from the total SNP
hitOWA <- OWAall.1[OWAall.1$SNP %in% hitOWA.1$SNP,]
str(hitOWA)
###get the table 
##combine the SNP and the Ind.SNP


alltable1 <- as.data.frame(table(OWAall.1[,c(10)]))
colnames(alltable1)[colnames(alltable1)=="Var1"] <- "Ind.SNP Combinations"
colnames(alltable1)[colnames(alltable1)=="Freq"] <- "SNP Totally detected"
hitOWAtable1 <- as.data.frame(table(hitOWA[,c(10)]))
colnames(hitOWAtable1)[colnames(hitOWAtable1)=="Var1"] <- "Ind.SNP Combinations"
colnames(hitOWAtable1)[colnames(hitOWAtable1)=="Freq"] <- "SNPs Hit With OWA Genes"
##merge the two table 
t.2 <- plyr::join_all(list(alltable1, hitOWAtable1), by="Ind.SNP Combinations")
###wite out the table 
write.csv(t.2,file="OWAtable/IndSNP.csv")


###imput the table and then making the bar images 

floall.1 <- read.csv("OWA2/OWA2.csv",row.names = 1)
floall.1 <- floall.1[!(floall.1$Ind.SNP=="122+35256im"),]
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

##get the hit with OWA 

hitflo.1 <- read.csv("OWA2/AHG.owa.hit.T.csv",row.names = 1)
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
colnames(hitflotable1)[colnames(hitflotable1)=="Freq"] <- "SNPs Hit With OWA Genes"
##merge the two table 
t.5 <- plyr::join_all(list(alltable1, hitflotable1), by="Non.missing")
str(t.5)
library(reshape2)
library(reshape)
t.16<- melt(t.5, id.vars = c("Non.missing"))
str(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With OWA Genes", "SNPs Totally Detected"))
write.csv(t.16, file="data/SNPs.M.csv")
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
t.16[t.16 == 0] <- NA
t.16$Non.missing <- factor(t.16$Non.missing,levels = c("O102+4322", "O112+3450", "O122+2646")) 
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
NON.mis <- ggplot(data=t.16, aes(x=Non.missing, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(labels=c("SNPs Hit with OWA Genes","SNPs Total Detected"))+
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
SNPcombs <- read.csv("OWAtable/IndSNP.combination.csv",header=T)
str(SNPcombs)

#colnames(t.15)[colnames(t.15)=="SNP.Totally.Detected.1"] <- "SNP.Totally.Detected"
#colnames(t.15)[colnames(t.15)=="SNPs.Hit.With.OWA.Genes.1"] <- "SNPs.Hit.With.OWA.Genes"
library(reshape2)
library(reshape)
t.18 <- melt(SNPcombs, id.vars = c("Type"))
str(t.18)
library(ggplot2)
t.18$variable<- factor(t.18$variable, levels= c("SNPs.Hit.With.OWA.Genes", "SNP.Totally.Detected"))

#t.18$Method <- factor(t.18$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
t.18[t.18 == 0] <- NA
t.18$Type <- factor(t.18$Type,levels = c("Imputed", "Non-missing", "Imputed and Non-missing"))
#t.18$Trait.name<- factor(t.18$Trait.name,levels = c("D", "W", "GW", "M"))

IM.NON.Over <- ggplot(data=t.18, aes(x=Type, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(labels=c("SNPs.Hit.With.OWA.Genes","SNPs.Totally.Detected"))+
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
floall.1 <- read.csv("OWA2/OWA2.csv",row.names = 1)
floall.1 <- floall.1[floall.1$Trait.name=="OWA19",]
floall.1 <- droplevels(floall.1)
str(floall.1)
levels
##select the method from mrMLM 
floall.1 <- floall.1[floall.1$Software=="mrMLMM",]
floall.1 <- droplevels(floall.1)
str(floall.1)
hitflo.1 <- read.csv("OWA2/AHG.owa.hit.T.csv")

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
#write.csv(t.3,file="Flotable/t.3.csv")
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
colnames(hitflotable1)[colnames(hitflotable1)=="Freq"] <- "SNPs Hit With OWA Genes"
##merge the two table 
t.19 <- plyr::join_all(list(alltable1, hitflotable1), by="PS")
t.19$Trait.name <- paste("OWA19")
###wite out the table 
write.csv(t.19,file="OWAtable/PS19.csv")

floall.1 <- read.csv("OWA2/OWA2.csv",row.names = 1)
floall.1 <- floall.1[floall.1$Trait.name=="OWA18",]
floall.1 <- droplevels(floall.1)
str(floall.1)
levels
##select the method from mrMLM 
floall.1 <- floall.1[floall.1$Software=="mrMLMM",]
floall.1 <- droplevels(floall.1)
str(floall.1)
hitflo.1 <- read.csv("OWA2/AHG.owa.hit.T.csv")

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

#write.csv(t.2,file="Flotable/t.2.csv")
t.3 <- t.1[t.1$SNP %in% t.2$SNP,]
t.5 <- t.1[!(t.1$SNP %in% t.2$SNP),]
#write.csv(t.3,file="Flotable/t.3.csv")
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
colnames(hitflotable1)[colnames(hitflotable1)=="Freq"] <- "SNPs Hit With OWA Genes"
##merge the two table 
t.18 <- plyr::join_all(list(alltable1, hitflotable1), by="PS")
t.18$Trait.name <- paste("OWA18")
###wite out the table 
write.csv(t.18,file="OWAtable/PS18.csv")

t.2 <- rbind(t.18,t.19)
###imput the table and then making the bar images 
###this table made by hnad "Ind.SNP.Combinations"
str(t.2)
library(reshape2)
library(reshape)
t.16<- melt(t.2, id.vars = c("PS","Trait.name"))
str(t.16)
#levels(t.16$PS)[levels(t.16$PS)=="PS_N;PS_Y"] <- "PS_N and PS_Y"
str(t.16)
library(ggplot2)
t.16$variable<- factor(t.16$variable, levels= c("SNPs Hit With OWA Genes", "SNP Totally Detected"))
#t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
#t.16[t.16 == 0] <- NA
t.16[t.16 == NA ] <- 0
t.16$PS <- factor(t.16$PS,levels = c("PS_N", "PS_Y", "PS_Y and PS_N"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))

PS <- ggplot(data=t.16, aes(x=t.16$PS, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="PS_N:No Population Structure\nPS_Y:Population Structure\nPS_Y and PS_N:Overlap",
                      labels=c("SNPs Hit With OWA Genes","SNPs.Totally.Detected"))+   
  geom_text(aes(label=t.16$value),size=4,position = position_stack(vjust = 0.5))+ facet_grid(~Trait.name)+
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
filename <- read.csv("OWA2/OWA2.csv",row.names = 1)
#levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="F116im+36088"] <- "F116+36088im"
levels(filename$Method)
str(filename)

library(tidyr)
library(dplyr)
##select the column need for study
t <- as.data.frame(ftable(filename[,c(12,15,16)]))
t.T.1 <- droplevels(t)
str(t.T.1)

###unit the SNP. method and trait.nmes 
t.T.1.u <-  unite_(t.T.1, "IndMT", c("Ind.SNP","Method", "Trait.name"),sep=";")
##change the column name Freq to  "SNPs.Total.detected"
colnames(t.T.1.u)[colnames(t.T.1.u)=="Freq"] <- "SNPs.Totally.detected"
###import the SNP hit with OWA genes 
hitOWA <- read.csv("OWA2/AHG.owa.hit.T.csv")
##select the SNPs hit with OWA but in the totole detected SNPs. 
hitflo <- filename[filename$SNP %in% hitOWA$SNP,]
hitflo <- droplevels(hitflo)
tflo <- as.data.frame(ftable(hitflo[,c(12,15,16)]))
nlevels(hitflo$SNP)
###unite the Method and Trait.name
###unit the SNP. method and trait.nmes in order to use to merge with the total detected SNPs
t.T.2.u <-  unite_(tflo, "IndMT", c("Ind.SNP", "Method", "Trait.name"),sep=";")
## change the column name Freq to "SNPs.hit.with.OWA.genes"
colnames(t.T.2.u)[colnames(t.T.2.u)=="Freq"] <- "SNPs.hit.with.OWA.genes"
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

t.15 <- t.14[t.14$Method=="MLM+SUPER" | t.14$Method=="CMLM+SUPER"|t.14$Method=="CMLM"|t.14$Method=="MLM"|t.14$Method=="GLM"|t.14$Method=="rrBLUP",]
#t.16 <- t.15[t.15$Ind.SNP=="F116+36088im",]
t.16 <- t.15
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.OWA.genes", "SNPs.Totally.detected"))
t.16$Ind.SNP <- factor(t.16$Ind.SNP,levels = c("122+35256im", "O102+4322", "O112+3450", "O122+2646"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
total.1 <- ggplot(data=t.16, aes(x=Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+ 
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.OWA.Genes","SNPs.Totally.Detected"))+   
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
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="Single-locus methods")

###Multiple methods 
t.15 <- t.14[t.14$Method=="FarmCPU" | t.14$Method=="mrMLM"|t.14$Method=="FASTmrMLM"|t.14$Method=="FASTmrEMMA"|t.14$Method=="pKWmEB"|
               t.14$Method=="pLARmEB"|t.14$Method=="ISIS EM-BLASSO",]
#t.16 <- t.15[t.15$Ind.SNP=="F116+36088im",]
t.16 <- t.15
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.OWA.genes", "SNPs.Totally.detected"))
t.16$Ind.SNP <- factor(t.16$Ind.SNP,levels = c("122+35256im", "O102+4322", "O112+3450", "O122+2646"))
#t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
total.1 <- ggplot(data=t.16, aes(x=Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.OWA.Genes","SNPs.Totally.Detected"))+   
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
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="Multiple-locus methods")


###getting the multiple methods 

floall.1 <- read.csv("OWA2/OWAuniquetable.csv")
floall.1 <- floall.1[floall.1$Trait.names=="OWA19",]
str(floall.1)
floall.1 <- droplevels(floall.1)
str(floall.1)
floallmethod1 <-  unite_(floall.1, "ThreeMeth", c("GAPIT","FarmCPU", "mrMLM"),sep=" ")
floallmethod.table1 <- as.data.frame(table(floallmethod1$ThreeMeth))
colnames(floallmethod.table1)[colnames(floallmethod.table1)=="Freq"] <- "SNPs.Totally.Detected"

hitflo.1 <- read.csv("OWA2/AHG.owa.hit.T.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo <- floall.1[floall.1$SNP %in% hitflo.1$SNP,]
str(hitflo)

hitflomethod <-  unite_(hitflo, "ThreeMeth", c("GAPIT","FarmCPU", "mrMLM"),sep=" ")
hitflomethod.table <- as.data.frame(table(hitflomethod$ThreeMeth))
colnames(hitflomethod.table)[colnames(hitflomethod.table)=="Freq"] <- "SNPs.Hit.With.OWA.Genes"

t.2 <- plyr::join_all(list(floallmethod.table1,hitflomethod.table), by="Var1")
t.2[is.na(t.2)] <- 0
write.csv(t.2, file="OWAtable/allmethod19.csv")
##separate MT to method and Trait.name after merge
#t.13 <- separate(t.2, "IndMT", c("Ind.SNP","Method", "Trait.name"),sep=";")
library(reshape2)
library(reshape)
t.14 <- melt(t.2, id.vars = c("Var1"))

str(t.14)
t.16 <- t.14
t.16[t.16 == 0] <- NA
t.16$variable<- factor(t.16$variable, levels= c("SNPs.Hit.With.OWA.Genes", "SNPs.Totally.Detected"))
allMeth <- ggplot(data=t.16, aes(x=reorder(Var1,-value), y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(labels=c("SNPs.Hit.With.OWA.Genes","SNPs.Totally.Detected"))+   
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
  labs(y="The number of QTNs",title="The number of SNPs from all methods")


###getting the multiple methods 

floall.1 <- read.csv("OWA2/OWAuniquetable.csv")
floall.2 <- floall.1[floall.1$Trait.names=="OWA18",]
str(floall.2)
str(floall.2)
floall.2 <- droplevels(floall.2)
str(floall.2)
floallmethod2 <-  unite_(floall.2, "ThreeMeth", c("GAPIT","FarmCPU", "mrMLM"),sep=" ")
floallmethod.table2 <- as.data.frame(table(floallmethod2$ThreeMeth))
colnames(floallmethod.table2)[colnames(floallmethod.table2)=="Freq"] <- "SNPs.Totally.Detected"

hitflo.1 <- read.csv("OWA2/AHG.owa.hit.T.csv")
###find the SNP hit with flowring genes from the total SNP
hitflo.3 <- floall.2[floall.2$SNP %in% hitflo.1$SNP,]
str(hitflo.3)

hitflomethod3 <-  unite_(hitflo.3, "ThreeMeth", c("GAPIT","FarmCPU", "mrMLM"),sep=" ")
hitflomethod.table3 <- as.data.frame(table(hitflomethod3$ThreeMeth))
colnames(hitflomethod.table3)[colnames(hitflomethod.table3)=="Freq"] <- "SNPs.Hit.With.OWA.Genes"

t.3 <- plyr::join_all(list(floallmethod.table2,hitflomethod.table3), by="Var1")
t.3[is.na(t.3)] <- 0
write.csv(t.3, file="OWAtable/allmethod18.csv")
##separate MT to method and Trait.name after merge
#t.13 <- separate(t.2, "IndMT", c("Ind.SNP","Method", "Trait.name"),sep=";")
library(reshape2)
library(reshape)
t.14 <- melt(t.3, id.vars = c("Var1"))

str(t.14)
t.16 <- t.14
t.16[t.16 == 0] <- NA
t.16$variable<- factor(t.16$variable, levels= c("SNPs.Hit.With.OWA.Genes", "SNPs.Totally.Detected"))
allMeth18 <- ggplot(data=t.16, aes(x=reorder(Var1,-value), y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(labels=c("SNPs.Hit.With.OWA.Genes","SNPs.Totally.Detected"))+   
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
  labs(y="The number of QTNs",title="The number of SNPs from all methods")


####getting all of the methods

