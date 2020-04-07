
library(ggplot2)

library(ggmap)
#install.packages("sp")
#install.packages("maps")
#install.packages("maptools")
library(sp)
library(maptools)
library(maps)

###install some packages 
##import the dataset 
library(devtools)
install_github('lchiffon/REmap')
library(REmap)
###this one for simple image of USA and 
tuskegee.mis <- read.csv("Charpter2/ggmaptuskegee.csv")
library(maps)
###USA
qplot(long,lat,data=tuskegee.mis,col="Red",size=8) + borders("state",size=0.5)
###Alabma 
qplot(long,lat,data=tuskegee.mis,col="Red",size=8) + borders("county","alabama",colour="grey70") + geom_point(colour = "Red",alpha=2)

####this one for ggmap like the image from the google map
####this one for ggmap like the image from the google map
install.packages("manipulate")
install.packages('maps')
install.packages("RgoogleMaps")
install.packages("maptools")
install.packages("sp")
install.packages("ggmap")
install.packages("ggplot2")
devtools::install_github("dkahle/ggmap", ref = "tidyup")

library(maps)
library(RgoogleMaps)
library(maptools)
library(sp)
library(ggmap)
library(ggplot2)
###
##USA <- get_map(location="state") if drawing the map of USA, it works with location ="state".

###if using ggmap, it need API key
### the way to get the API, it to "?get_googlemap" to get the web, and then regist, it will get 300$ for 12 months free API key. 
###https://developers.google.com/maps/documentation/maps-static/intro, 
###https://developers.google.com/maps/documentation/maps-static/dev-guide,
###https://developers.google.com/maps/documentation/maps-static/get-api-key, 
###https://developers.google.com/maps/documentation/maps-static/usage-and-billing ggmap, register_google
###and the go the API key and then paste it in the key=""

register_google(key = "AIzaSyB1v5A-0Qg2nD90Rv16rOqqG2b_z2KvoBY")
has_google_key()
###making the image ready 

install.packages("readr")
library(readr)

install.packages("proj4")
library(proj4)

install.packages("magick")
library(magick)

#install.packages("ggmap")
library(ggmap)

#install.packages("ggimage")
library(ggimage)  
mis <-image_transparent(image_read("Miscan1.png"), 'white')
#mis <- image_colorize(mis, 100, "darkgreen")
image_write(mis, path = "Miscan2.png", format = "png")
###import the dataset 

add <- read.csv("ggmaptuskegee.csv", header = TRUE, sep = ",", dec = ".")
add2 <-add %>% 
  mutate(Image = case_when(city == "Tuskegee" ~ "Miscan2.png"))
###this one is to get the image of the alabama

mapt <- get_googlemap(center = c(-88, 32), zoom = 6, source = "google", maptype = 'terrain')
ggmap(mapt) +
  geom_image(aes(x = long, y = lat, image=Image), data = add2, size = 0.03)

###this one for the whole America
mapt <- get_googlemap(center = c(-97, 39), zoom = 4, source = "google", maptype = 'terrain')

ggmap(mapt) +
  geom_image(aes(x = long, y = lat, image=Image), data = add2, size = 0.03)

###this one using the mapview and like the google map, which can change the size of the image. 
library(ggmap)
library(sp)
library(mapview)
install.packages("mapview") ###restart R section, it works
install.packages("ggimage")

#add <- as.data.frame(read.csv("hainan.csv", header = TRUE, sep = ",", dec = "."))
coordinates(add) <- ~ long + lat
#plot(add)
proj4string(add) <- CRS("+init=epsg:4326")
###if does not work in the macBook pro, it will need download "XQuartz" 
mapview(add)
