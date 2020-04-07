
###GLM
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FD_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FD_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM.W <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FW_1.GWAS.Results.csv")
str(GLM.W)
GLM.W<- GLM.W[1:4]
colnames(GLM.W)[colnames(GLM.W)=="P.value"] <- "GLM.FW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM.GW <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GFW_1.GWAS.Results.csv")
str(GLM.GW)
GLM.GW<- GLM.GW[1:4]
colnames(GLM.GW)[colnames(GLM.GW)=="P.value"] <- "GLM.GFW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM.M <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FM_1.GWAS.Results.csv")
str(GLM.M)
GLM.M<- GLM.M[1:4]
colnames(GLM.M)[colnames(GLM.M)=="P.value"] <- "GLM.FM_1"
plot <-  plyr::join_all(list(GLM,GLM.W,GLM.GW,GLM.M), by="SNP")
str(plot)
plot <- plot[,c(1:4,7,10,13)]
str(plot)
###plot QQ plot inone image
setwd("GLMCVFD1image")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=8.85e-7,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(8.85e-7,9.83e-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(8.85e-7,9.83e-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)



###FarmCPU
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("GAPIT005/FarmCPU0.05/FarmCPU.FD_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FD_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU.W <- read.csv("GAPIT005/FarmCPU0.05/FarmCPU.FW_1.GWAS.Results.csv")
str(FarmCPU.W)
FarmCPU.W<- FarmCPU.W[1:4]
colnames(FarmCPU.W)[colnames(FarmCPU.W)=="P.value"] <- "FarmCPU.FW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU.GW <- read.csv("GAPIT005/FarmCPU0.05/FarmCPU.GFW_1.GWAS.Results.csv")
str(FarmCPU.GW)
FarmCPU.GW<- FarmCPU.GW[1:4]
colnames(FarmCPU.GW)[colnames(FarmCPU.GW)=="P.value"] <- "FarmCPU.GFW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU.M <- read.csv("GAPIT005/FarmCPU0.05/FarmCPU.FM_1.GWAS.Results.csv")
str(FarmCPU.M)
FarmCPU.M<- FarmCPU.M[1:4]
colnames(FarmCPU.M)[colnames(FarmCPU.M)=="P.value"] <- "FarmCPU.FM_1"
plot <-  plyr::join_all(list(FarmCPU,FarmCPU.W,FarmCPU.GW,FarmCPU.M), by="Name")
colnames(plot)[colnames(plot)=="Name"] <- "SNP"
str(plot)
plot <- plot[,c(1:4,7,10,13)]
str(plot)
###plot QQ plot inone image
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("FarmCPUCVFD1image1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
F1 <- CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=8.85e-7,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
F2 <- CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(8.85e-7,9.83e-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
F3 <- CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(8.85e-7,9.83e-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)




###FarmCPU
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("GAPIT005/FarmCPUCV0.05/FarmCPU.FD_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU_PS.FD_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU.W <- read.csv("GAPIT005/FarmCPUCV0.05/FarmCPU.FW_1.GWAS.Results.csv")
str(FarmCPU.W)
FarmCPU.W<- FarmCPU.W[1:4]
colnames(FarmCPU.W)[colnames(FarmCPU.W)=="P.value"] <- "FarmCPU_PS.FW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU.GW <- read.csv("GAPIT005/FarmCPUCV0.05/FarmCPU.GFW_1.GWAS.Results.csv")
str(FarmCPU.GW)
FarmCPU.GW<- FarmCPU.GW[1:4]
colnames(FarmCPU.GW)[colnames(FarmCPU.GW)=="P.value"] <- "FarmCPU_PS.GFW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU.M <- read.csv("GAPIT005/FarmCPUCV0.05/FarmCPU.FM_1.GWAS.Results.csv")
str(FarmCPU.M)
FarmCPU.M<- FarmCPU.M[1:4]
colnames(FarmCPU.M)[colnames(FarmCPU.M)=="P.value"] <- "FarmCPU_PS.FM_1"
plot <-  plyr::join_all(list(FarmCPU,FarmCPU.W,FarmCPU.GW,FarmCPU.M), by="Name")
str(plot)
plot <- plot[,c(1:4,7,10,13)]
str(plot)
###plot QQ plot inone image
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("FarmCPUCVFD1image1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
F4 <- CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=8.85e-7,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
F5 <- CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(8.85e-7,9.83e-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
F6 <- CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(8.85e-7,9.83e-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

FCVNOCV <- ggarrange(F1,F4, labels = c("A", "B"), ncol = 2, nrow = 1,common.legend = TRUE, legend="bottom" )



###MLMSUPER
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
MLMSUPER <- read.csv("GAPIT005/MLMSUPPER0.05/GAPIT.SUPER.FD_1.GWAS.Results.csv")
MLMSUPER <- MLMSUPER[1:4]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "MLMSUPER.FD_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
MLMSUPER.W <- read.csv("GAPIT005/MLMSUPPER0.05/GAPIT.SUPER.FW_1.GWAS.Results.csv")
str(MLMSUPER.W)
MLMSUPER.W<- MLMSUPER.W[1:4]
colnames(MLMSUPER.W)[colnames(MLMSUPER.W)=="P.value"] <- "MLMSUPER.FW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
MLMSUPER.GW <- read.csv("GAPIT005/MLMSUPPER0.05/GAPIT.SUPER.GFW_1.GWAS.Results.csv")
str(MLMSUPER.GW)
MLMSUPER.GW<- MLMSUPER.GW[1:4]
colnames(MLMSUPER.GW)[colnames(MLMSUPER.GW)=="P.value"] <- "MLMSUPER.GFW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
MLMSUPER.M <- read.csv("GAPIT005/MLMSUPPER0.05/GAPIT.SUPER.FM_1.GWAS.Results.csv")
str(MLMSUPER.M)
MLMSUPER.M<- MLMSUPER.M[1:4]
colnames(MLMSUPER.M)[colnames(MLMSUPER.M)=="P.value"] <- "MLMSUPER.FM_1"
plot <-  plyr::join_all(list(MLMSUPER,MLMSUPER.W,MLMSUPER.GW,MLMSUPER.M), by="SNP")
str(plot)
plot <- plot[,c(1:4,7,10,13)]
str(plot)
###plot QQ plot inone image
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("MLMSUPERCVFD1image1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=8.85e-7,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(8.85e-7,9.83e-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(8.85e-7,9.83e-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###CMLMSUPER
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
CMLMSUPER <- read.csv("GAPIT005/CMLMSUPPER0.05/GAPIT.SUPER.FD_1.GWAS.Results.csv")
CMLMSUPER <- CMLMSUPER[1:4]
colnames(CMLMSUPER)[colnames(CMLMSUPER)=="P.value"] <- "CMLMSUPER.FD_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
CMLMSUPER.W <- read.csv("GAPIT005/CMLMSUPPER0.05/GAPIT.SUPER.FW_1.GWAS.Results.csv")
str(CMLMSUPER.W)
CMLMSUPER.W<- CMLMSUPER.W[1:4]
colnames(CMLMSUPER.W)[colnames(CMLMSUPER.W)=="P.value"] <- "CMLMSUPER.FW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
CMLMSUPER.GW <- read.csv("GAPIT005/CMLMSUPPER0.05/GAPIT.SUPER.GFW_1.GWAS.Results.csv")
str(CMLMSUPER.GW)
CMLMSUPER.GW<- CMLMSUPER.GW[1:4]
colnames(CMLMSUPER.GW)[colnames(CMLMSUPER.GW)=="P.value"] <- "CMLMSUPER.GFW_1"

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
CMLMSUPER.M <- read.csv("GAPIT005/CMLMSUPPER0.05/GAPIT.SUPER.FM_1.GWAS.Results.csv")
str(CMLMSUPER.M)
CMLMSUPER.M<- CMLMSUPER.M[1:4]
colnames(CMLMSUPER.M)[colnames(CMLMSUPER.M)=="P.value"] <- "CMLMSUPER.FM_1"
plot <-  plyr::join_all(list(CMLMSUPER,CMLMSUPER.W,CMLMSUPER.GW,CMLMSUPER.M), by="SNP")
str(plot)
plot <- plot[,c(1:4,7,10,13)]
str(plot)
###plot QQ plot inone image
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("CMLMSUPERCVFD1image1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=8.85e-7,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(8.85e-7,9.83e-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(8.85e-7,9.83e-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
