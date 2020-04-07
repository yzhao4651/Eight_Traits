###this one used for search the genes around SNPs

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")

OWA2 <- read.csv("OWA2/OWAunique.csv")
str(OWA2)
levels(OWA2$SNP)

Chr11.gff3.txt

install.packages("BiocManager")
BiocManager::install("rtracklayer")
BiocManager::install("GenomicFeatures")
BiocManager::install("IRanges")
BiocManager::install("AnnotationDbi")
BiocManager::install("GenomicFeatures")
BiocManager::install("AnnotationDbi")
source("https://bioconductor.org/biocLite.R")
BiocManager::install("biocLite")
BiocManager::install("Repitools")
library(Repitools)
library(AnnotationDbi)
library(rtracklayer)
library(GenomicFeatures)
library(IRanges)
txdb <- makeTxDbFromGFF("Sbicolor_454_v3.1.1.gene.gff3",
                        format = "gff3",
                        dataSource = "Phytosome 12",
                        organism = "Sorghum bicolor")

mySNPs <- OWA2
source("Function/Findneargenes.R")
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    saved_genes[[i]] <- mygenes$tx_name
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        Gene = unlist(saved_genes))
  return(gene_df)
}

gene_df.genes <- Findneargenes(OWA2, 1e4)
str(gene_df.genes)
nlevels(gene_df.genes$Gene)
#seqinfo(txdb)
#To save the txdb database for later and avoid having to recreate it every time we use it, we can use saveDb() and, later, loadDb()
#BiocManager::install("AnnotationDbi")
#library(AnnotationDbi)
#saveDb(txdb, 'txdb.Sbicolor_454_v3.1.1')
#Once the database is stored to disk, it can be reloaded later, almost instantly.
#library(AnnotationDbi)
#txdb = loadDb(file = 'txdb.Sbicolor_454_v3.1.1')
#genes(txdb)
#seqinfo(txdb)
###install.packags library(IRanges)
#BiocManager::install("IRanges")
### fucntion to get the neargenes 
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
### this one for flowering traits
#full_data.flo <- read.csv("data/flo-106.116.122.106im.csv")
#str(full_data.flo)
#levels(full_data.flo$SNP)
#nlevels(full_data.flo$SNP)
##this function need input dataset and search_radius 
#source("https://bioconductor.org/biocLite.R")
#biocLite("Repitools")
library(Repitools)
columns(txdb)
mySNPs <- OWA2
source("Function/Findneargenes.R")
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
      t <- annoGR2DF(mygenes) 
    saved_genes[[i]] <- t$chr
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        chr = unlist(saved_genes))
  return(gene_df)
}


gene_df.chr <- Findneargenes(OWA2, 1e4)
str(gene_df.chr)

Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)

    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$start
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        start = unlist(saved_genes))
  return(gene_df)
}

gene_df.start <- Findneargenes(OWA2, 1e4)
str(gene_df.start)


Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$end
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        end = unlist(saved_genes))
  return(gene_df)
}

gene_df.end <- Findneargenes(OWA2, 1e4)
str(gene_df.end)
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$width
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        width = unlist(saved_genes))
  return(gene_df)
}

gene_df.width <- Findneargenes(OWA2, 1e4)
str(gene_df.width)

Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$strand
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        strand = unlist(saved_genes))
  return(gene_df)
}

gene_df.strand <- Findneargenes(OWA2, 1e4)
str(gene_df.strand)

###combine all of them together 
gene_df10.17 <- do.call("cbind",list(gene_df.genes,gene_df.chr, gene_df.start, gene_df.end, gene_df.width, gene_df.strand))
str(gene_df10.17)

gene_df10.17.a <- gene_df10.17[,c(1,2,4,6,8,10,12)]
gene_df10.17.a <- droplevels(gene_df10.17.a)

nlevels(gene_df10.17.a$SNP) ##111
nlevels(gene_df10.17.a$Gene)###417

write.csv(gene_df10.17.a,file="OWA2/AHG.OWA2.T.hitandnohit.csv")

##merge all of data withunique

OWA2.T <- merge(OWA2, gene_df10.17.a,by="SNP",all.x=TRUE)
OWA2.T <- droplevels(OWA2.T)
nlevels(OWA2.T$SNP)##111
nlevels(OWA2.T$Gene)##417
nlevels(OWA2.T$SNP)###135
###get all of the no any hit SNP
OWA2.T.nohhit <- OWA2.T[is.na(OWA2.T$Gene),]
OWA2.T.nohhit <- droplevels(OWA2.T.nohhit)
nlevels(OWA2.T.nohhit$SNP)###24
write.csv(OWA2.T.nohhit,file="OWA2/AHG.OWA2.T.nohit.csv")

###all of them with hit gene
OWA2.T.hit <- OWA2.T[!(is.na(OWA2.T$Gene)),]
str(OWA2.T.hit)
OWA2.T.hit$dis_start.positon <- OWA2.T.hit$start - OWA2.T.hit$Position
OWA2.T.hit$dis_positon.end <- OWA2.T.hit$Position- OWA2.T.hit$end
OWA2.T.hit <- droplevels(OWA2.T.hit)
nlevels(OWA2.T.hit$SNP) ###111
nlevels(OWA2.T.hit$Gene)###417

write.csv(OWA2.T.hit,file="OWA2/AHG.OWA2.T.hit.csv")
###and then get gene anotation 
annotations <- read.delim("Sbicolor_454_v3_1_1_annotation_info.txt",
                          stringsAsFactors = FALSE)
head(annotations)
str(annotations)
myrows <- match(OWA2.T.hit$Gene, annotations$transcriptName)
allhitsgenes.owa <- annotations[myrows,]
head(allhitsgenes.owa)
###combine two dataset in order to get the SNP with gene anotation in one file. 
AHG.owa.1 <- cbind(OWA2.T.hit,allhitsgenes.owa)
colnames(AHG.owa.1)[colnames(AHG.owa.1)=="Gene"] <- "Gene"

library(dplyr)
install.packages("dplyr")
install.packages("tidyverse")
library(tidyverse)
##seperate the Gene with number and then delete the duplicated genes
AHG.owa.2 <- separate(AHG.owa.1, Gene, c("Gene1","Gene2","Gene.num"), "\\.") 
AHG.owa.3 <- unite(AHG.owa.2,"SNP.Gene",SNP,Gene1,Gene2,sep="/")
AHG.owa.4<- AHG.owa.3[!(duplicated(AHG.owa.3$SNP.Gene)),]
AHG.owa.5 <- separate(AHG.owa.4, SNP.Gene, c("SNP","Gene1","Gene2"), "/")
AHG.owa.5 <- unite(AHG.owa.5, Gene, c("Gene1","Gene2"), sep=".")
AHG.owa.5 <- droplevels(AHG.owa.5)
nlevels(as.factor(AHG.owa.5$SNP)) ###111
nlevels(as.factor(AHG.owa.5$Gene))###317

write.csv(AHG.owa.5,file="OWA2/AHG.OWA2.anota.csv")
###select the Genes without studying: which means the X.pac is NA 
owa.hit.gene <- read.csv("OWA2/AHG.OWA2.anota.csv")
str(owa.hit.gene)
nlevels(as.factor(owa.hit.gene$Gene))

###this one works even this is hit with sobic genes, but many of them has no any study 80 has no any function found 
owa.NOX.pacId <- owa.hit.gene[is.na(owa.hit.gene$X.pacId),]
owa.NOX.pacId <- droplevels(owa.NOX.pacId)
nlevels(owa.NOX.pacId$SNP)###10
nlevels(owa.NOX.pacId$Gene)###30
write.csv(owa.NOX.pacId,file="OWA2/AHG.owa.NOX.pacId.csv")

## with some hit 
owa.X.pacId <- owa.hit.gene[!is.na(owa.hit.gene$X.pacId),]
owa.X.pacId <- droplevels(owa.X.pacId)
nlevels(owa.X.pacId$SNP)###101
nlevels(owa.X.pacId$Gene)###287
write.csv(owa.X.pacId,file="OWA2/AHG.owa.X.pacId.csv")

###check how many hit with ara 
owa.NO.ara <- owa.X.pacId[owa.X.pacId$Best.hit.arabi.name=="",]
owa.NO.ara <- droplevels(owa.NO.ara)
names(owa.NO.ara)
nlevels(owa.NO.ara$SNP)###25
nlevels(owa.NO.ara$Gene)###35
write.csv(owa.NO.ara,file="OWA2/AHG.owa.X.pacId.NO.ara.csv")


###check how many hit with ara 
owa.X.pacId.ara <- owa.X.pacId[!(owa.X.pacId$Best.hit.arabi.name==""),]
owa.X.pacId.ara <- droplevels(owa.X.pacId.ara)
names(owa.X.pacId.ara)
nlevels(owa.X.pacId.ara$SNP)###100
nlevels(owa.X.pacId.ara$Gene)###252
write.csv(owa.X.pacId.ara,file="OWA2/AHG.owa.X.pacId.ara.csv")
###check the GO information
### and then using the data file to find the SNPs hit GO term 
owa.X.pacId.GO <- owa.X.pacId[!(owa.X.pacId$GO==""),]
owa.X.pacId.GO <- droplevels(owa.X.pacId.GO)
names(owa.X.pacId.GO)
nlevels(owa.X.pacId.GO$SNP)###89
nlevels(owa.X.pacId.GO$Gene)###162
write.csv(owa.X.pacId.GO,file="OWA2/AHG.owa.X.pacI.GO.csv")
#### rice 
owa.X.pacId.rice <- owa.X.pacId[!(owa.X.pacId$Best.hit.rice.name==""),]
owa.X.pacId.rice <- droplevels(owa.X.pacId.rice)
names(owa.X.pacId.rice)
nlevels(owa.X.pacId.rice$SNP)###101
nlevels(owa.X.pacId.rice$Gene)###262
write.csv(owa.X.pacId.rice,file="OWA2/AHG.owa.X.pacI.rice.csv")

###load the sorghum data base 
#this library can not be used in MacBook Pro
#install.packages("RODBC")
#library(RODBC)
#channel <- odbcConnectAccess2007("Sorghum_Cold_DB")

##using this one to check if there is cold genes in three species: sorgum
###this one for searching the flowering time genes 
###sepate the gene with .1 using function 
library(tidyr)
AHG.owa.2 <- separate(AHG.owa.1, Gene, c("Gene1","Gene2","Gene.num"), "\\.") 
### and then remove the same genes and same SNP
AHG.owa.3 <- AHG.owa.2[!duplicated(AHG.owa.2[,c("SNP","Gene2")]),]
AHG.owa.3 <- unite(AHG.owa.3,"Gene",Gene1,Gene2,sep=".")
str(AHG.owa.3)
AHG.owa.4 <- separate(AHG.owa.3, Best.hit.arabi.name, c("Best.hit.arabi.name","Best.hit.arabi.name.num"), "\\.") 
###import the candidate gene from ara
Ara.gene.accli <- read.csv("OWA/ara.cold.accli.csv", sep='\t')
Ara.gene.resp.c <- read.csv("OWA/ara.respon.cold.csv", sep='\t')

str(Ara.gene.resp.c)
str(Ara.gene.accli)

##hit with owawering time genes with ara

AHG.owa.1.ara.hit <- AHG.owa.1[AHG.owa.1$Best.hit.arabi.name %in% Ara.gene.accli$Gene.Model,]
AHG.owa.1.ara.hit.c <- AHG.owa.1[AHG.owa.1$Best.hit.arabi.name %in% Ara.gene.resp.c$Gene.Model,]

AHG.owa.1.ara.hit.c$hit.gene <- paste("arabidopsis")
names(AHG.owa.1.ara.hit.c)
###hit with cold with rice 
library(tidyr)
AHG.owa.2 <- separate(AHG.owa.1, Gene, c("Gene1","Gene2","Gene.num"), "\\.") 
AHG.owa.3 <- AHG.owa.2[!duplicated(AHG.owa.2[,c("SNP","Gene2")]),]
AHG.owa.3 <- unite(AHG.owa.3,"Gene",Gene1,Gene2,sep=".")
AHG.owa.5 <- separate(AHG.owa.3, Best.hit.rice.name, c("Best.hit.rice.name","Best.hit.rice.name.num"), "\\.") 
str(AHG.owa.5)
rice.gene <- read.csv("OWA/OWA.rice.csv")
str(rice.gene)
AHG.owa.1.Rice.hit <- AHG.owa.5[AHG.owa.5$Best.hit.rice.name %in% rice.gene$genes,]
AHG.owa.1.Rice.hit$hit.gene <- paste("oryza(rice)")
str(AHG.owa.1.Rice.hit)
AHG.owa.1.Rice.hit <- AHG.owa.1.Rice.hit[,-c(15,37)]
names(AHG.owa.1.Rice.hit)

###hit with cold with sorghum
##first need change the gene name to the name used in this project 
##import the data need change the name
Surghum.name.need.change <- read.csv("OWA/OWA.sorghum.go.csv")
Surghum.bothnamefile <- read.csv("OWA/SB.GI.exchange.csv")
Surghum.bothnamefile <- Surghum.bothnamefile[1:2]
##get the file with right name
Sorghum.right.name.1 <- Surghum.bothnamefile[Surghum.bothnamefile$SB %in% Surghum.name.need.change$Gene.1,]
Sorghum.right.name.2 <- Surghum.bothnamefile[Surghum.bothnamefile$SB %in% Surghum.name.need.change$gene.2,]
Sorghum.right.name.3 <- Surghum.bothnamefile[Surghum.bothnamefile$SB %in% Surghum.name.need.change$gene.3,]
Sorghum.right.name.4 <- Surghum.bothnamefile[Surghum.bothnamefile$SB %in% Surghum.name.need.change$gene.4,]
Sorghum.right.name.5 <- Surghum.bothnamefile[Surghum.bothnamefile$SB %in% Surghum.name.need.change$gene.5,]
##write out the result
write.csv(Sorghum.right.name.1,file="OWA2/Sorghum.right.name.1.csv")
write.csv(Sorghum.right.name.2,file="OWA2/Sorghum.right.name.2.csv")
write.csv(Sorghum.right.name.3,file="OWA2/Sorghum.right.name.3.csv")
write.csv(Sorghum.right.name.4,file="OWA2/Sorghum.right.name.4.csv")
write.csv(Sorghum.right.name.5,file="OWA2/Sorghum.right.name.5.csv")
##this one use for merge with different column names and different number of rows
#Sorghum.1.2 <- do.call("merge", c(lapply(list(Sorghum.right.name.1, Sorghum.right.name.2), data.frame, row.names=NULL), by = 0, all = TRUE))[-1]
#Sorghum.1.2.3 <- do.call("merge", c(lapply(list(Sorghum.1.2, Sorghum.right.name.3), data.frame, row.names=NULL), by = 0, all = TRUE))[-1]
#Sorghum.1.2.3.4<- do.call("merge", c(lapply(list(Sorghum.1.2.3, Sorghum.right.name.4), data.frame, row.names=NULL), by = 0, all = TRUE))[-1]
#Sorghum.all <- do.call("merge", c(lapply(list(Sorghum.1.2.3.4, Sorghum.right.name.5), data.frame, row.names=NULL), by = 0, all = TRUE))[-1]
#Sorghum.1.2 <- merge(data.frame(Sorghum.right.name.1, row.names=NULL), data.frame(Sorghum.right.name.2, row.names=NULL), by = 0, all = TRUE)[-1]

###report the sorghum genes 
Sorghum.genes <- read.csv("OWA/Sorghum.genes.cold.csv")
str(Sorghum.genes)
levels(Sorghum.genes$gene.ID)

AHG.owa.2 <- separate(AHG.owa.1, Gene, c("Gene1","Gene2","Gene.num"), "\\.") 
AHG.owa.3 <- AHG.owa.2[!duplicated(AHG.owa.2[,c("SNP","Gene2")]),]
AHG.owa.3 <- unite(AHG.owa.3,"Gene",Gene1,Gene2,sep=".")
str(AHG.owa.3)
AHG.owa.1.Sorghum.hit <- AHG.owa.3[AHG.owa.3$Gene %in% Sorghum.genes$gene.ID,]
AHG.owa.1.Sorghum.hit$hit.gene <- paste("Sorghum")
AHG.owa.1.Sorghum.hit <- AHG.owa.1.Sorghum.hit[,-c(15)]
str(AHG.owa.1.Sorghum.hit)
names(AHG.owa.1.Sorghum.hit)

###find if there are some SNP hit with Miscanthus 

MisOWAHongxu <- read.csv("OWA/MisOWASNP.csv")
str(MisOWAHongxu)
MisOWAHongxu$V1 <- make.names(MisOWAHongxu$V1)
OWA2 <- read.csv("OWA2/OWAunique.csv")
str(OWA2)
### 
AHG.owa.1.Mis.hit.1 <- AHG.owa.1[MisOWAHongxu$V1 %in% AHG.owa.1$SNP,]
AHG.owa.1.Mis.hit.2 <- OWA2[MisOWAHongxu$V1 %in% OWA2$SNP,]
AHG.owa.1.Mis.hit.1 <- droplevels(AHG.owa.1.Mis.hit.1)
AHG.owa.1.Mis.hit.2 <- droplevels(AHG.owa.1.Mis.hit.2)

####search the Sorghum genes from the Miscanthus hit 
library(tidyr)
library(readxl)
allsnp <- read_excel("data/allsnp.xlsx")
allsnp.1 <- data.frame(allsnp[,c(1:2,9:10)])
str(allsnp.1)
Mis.sorghum <- read.csv("OWA/OWA.miscanthus.csv")
str(Mis.sorghum)
Mis.sorghum.1 <- allsnp.1[allsnp.1$Marker.name  %in%  Mis.sorghum$Miscanthus.SNPe,]

###merge Mssorghum.1 with Mis.sorghum 
Mis.sorghum.all <- merge(Mis.sorghum,Mis.sorghum.1, by.x="Miscanthus.SNPe",by.y="Marker.name",x.all=TRUE)
str(Mis.sorghum.all)
###write out the data set with the position and location 

write.csv(Mis.sorghum.all, file="OWA2/OWAHONGXU.csv")
###check if there are some names that same to us 
OWA2 <- read.csv("OWA2/OWAunique.csv")
str(OWA2)
###there is no any one is hit with the OWA from Hongxu Dong
AHG.owa.1.Mis.hit <- OWA2[OWA2$SNP %in% Mis.sorghum.all$Original.marker.name, ]

###check if there is one hit with the same arabidropsis 
library(tidyr)
AHG.owa.2 <- separate(AHG.owa.1, Gene, c("Gene1","Gene2","Gene.num"), "\\.") 
### and then remove the same genes and same SNP
AHG.owa.3 <- AHG.owa.2[!duplicated(AHG.owa.2[,c("SNP","Gene2")]),]
AHG.owa.3 <- unite(AHG.owa.3,"Gene",Gene1,Gene2,sep=".")
str(AHG.owa.3)
AHG.owa.1.Mis.Sorghum.hit <- AHG.owa.3[AHG.owa.3$Gene %in% Mis.sorghum$Nearest.S..bicolor.gene,]
AHG.owa.1.Mis.Sorghum.hit$hit.gene <- paste("Mis.Sorghum")
AHG.owa.1.Mis.Sorghum.hit <- AHG.owa.1.Mis.Sorghum.hit[,-c(15)]

str(AHG.owa.1.Mis.Sorghum.hit)
names(AHG.owa.1.Mis.Sorghum.hit)

####search the arabidopisis genes from the Miscanthus hit 
AHG.owa.4 <- separate(AHG.owa.3, Best.hit.arabi.name, c("Best.hit.arabi.name","Best.hit.arabi.name.num"), "\\.") 
AHG.owa.4 <- AHG.owa.4[AHG.owa.4$Best.hit.arabi.name !="",]
Mis.sorghum <- read.csv("OWA/OWA.miscanthus.csv")
str(Mis.sorghum)
AHG.owa.1.Mis.ara.hit <- AHG.owa.4[AHG.owa.4$Best.hit.arabi.name %in% Mis.sorghum$Arabidopsis.gene,]
AHG.owa.1.Mis.ara.hit$hit.gene <- paste("Mis.ara")
names(AHG.owa.1.Mis.ara.hit)
AHG.owa.1.Mis.ara.hit <- AHG.owa.1.Mis.ara.hit[,-c(15,34)]

str(AHG.owa.1.Mis.ara.hit)
names(AHG.owa.1.Mis.ara.hit)
###search the rice species 
#Mis.sorghum <- read.csv("OWA/OWA.miscanthus.csv")
#str(Mis.sorghum)
#AHG.owa.1.Mis.Sorghum.hit <- AHG.owa.5[AHG.owa.5$Best.hit.rice.name %in% Mis.sorghum$Nearest.S..bicolor.gene,]
#AHG.owa.1.Mis.Sorghum.hit$hit.gene <- paste("Mis.rice")
                                   
###also can use this one "plyr::rbind.fill(mzList)"
AHG.owa.hit.T <- do.call("rbind",list(AHG.owa.1.Rice.hit, AHG.owa.1.Sorghum.hit, AHG.owa.1.Mis.Sorghum.hit, AHG.owa.1.Mis.ara.hit,
                                      AHG.owa.1.ara.hit.c))
all(names(AHG.owa.1.Rice.hit) %in% names(AHG.owa.1.Sorghum.hit))
all(names(AHG.owa.1.Rice.hit) %in% names(AHG.owa.1.Sorghum.hit))

all(names(AHG.owa.1.Rice.hit) %in% names(AHG.owa.1.Mis.Sorghum.hit))
all(names(AHG.owa.1.Rice.hit) %in% names(AHG.owa.1.Mis.ara.hit))
all(names(AHG.owa.1.Rice.hit) %in% names(AHG.owa.1.ara.hit.c))




write.csv(AHG.owa.hit.T, file="OWA2/AHG.owa.hit.T.csv")
###find the Sorghum genes 
Sorghum.genes <- read.csv("OWA/Sorghum.genes.cold.csv")

AHG.owa.hit.T.Sorghum <- Sorghum.genes[Sorghum.genes$gene.ID %in% AHG.owa.hit.T$Gene,]

write.csv(AHG.owa.hit.T.Sorghum,file="OWA2/AHG.owa.hit.T.Sorghum.csv")




