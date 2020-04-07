##import the dataset 

all.chapter2 <- read.csv("Chapter2UP/Supporting Information Table S1.Ch2.csv")

need.yld <- read.csv("OWA2/used accession to find their yld.csv")

all.chapter2.need <- all.chapter2[all.chapter2$Accession.. %in% need.yld$Accession,]
write.csv(all.chapter2.need,file="OWA2/OWA2.yeld.csv")