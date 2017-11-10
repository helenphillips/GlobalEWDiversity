
library(dplyr)
library(maptools)
library(maps)
library(plyr)
library(dplyr)


########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

##########################################################
## Load in data
##########################################################

data_in <-"0_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
loadin <- files[grep(date, files)]
loadinspecies <- loadin[grep("species_", loadin)]

species <- read.csv(file.path(data_in, loadinspecies))



data_in <-"3_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
loadin <- files[grep(date, files)]

sites <- read.csv(file.path(data_in, loadin))



############
## Where is the data
############


coord<-aggregate(cbind(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees), list(sites$Study_Name), mean)
## Five don't have coordinates yet
coord <- coord[complete.cases(coord),]
coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")
dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")


## Species data
coord2<-aggregate(cbind(species$Longitude__Decimal_Degrees, species$Latitude__decimal_degrees), list(species$Study_Name), mean)
## Five don't have coordinates yet
coord2 <- coord2[complete.cases(coord2),]
coord2$X<-coord2$Group.1
coord2<-coord2[2:4]
names(coord2)<-c("Long", "Lat", "X")
dsSPDF2<-SpatialPointsDataFrame(coord2[,1:2], data.frame(coord2[,1:3]))
proj4string(dsSPDF2)<-CRS("+proj=longlat")


pdf(file = file.path(figures, "Map_alldata_showingSpeciesLevel.pdf"), height = 4)
#jpeg(filename = file.path(figures, "Map_alldata.jpg"), quality = 100, res = 300, width = 2000, height = 2000)
mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)
points(dsSPDF2, col="red", bg="red", cex= 1, pch=19)
dev.off()



###### Invasives
tapply(species$Study_Name, species$Native.Nonnative, function(x) length(unique(x)))


natives <- species[species$Native.Nonnative == "Native",]
coord3<-aggregate(cbind(natives$Longitude__Decimal_Degrees, natives$Latitude__decimal_degrees), list(natives$Study_Name), mean)
## Five don't have coordinates yet
coord3 <- coord3[complete.cases(coord3),]
coord3$X<-coord3$Group.1
coord3<-coord3[2:4]
names(coord3)<-c("Long", "Lat", "X")
dsSPDF3<-SpatialPointsDataFrame(coord3[,1:2], data.frame(coord3[,1:3]))
proj4string(dsSPDF3)<-CRS("+proj=longlat")

nonnatives <- species[species$Native.Nonnative == "Non-native",]
coord4<-aggregate(cbind(nonnatives$Longitude__Decimal_Degrees, nonnatives$Latitude__decimal_degrees), list(nonnatives$Study_Name), mean)
## Five don't have coordinates yet
coord4 <- coord4[complete.cases(coord4),]
coord4$X<-coord4$Group.1
coord4<-coord4[2:4]
names(coord4)<-c("Long", "Lat", "X")
dsSPDF4<-SpatialPointsDataFrame(coord4[,1:2], data.frame(coord4[,1:3]))
proj4string(dsSPDF4)<-CRS("+proj=longlat")


pdf(file = file.path(figures, "Map_nativesnonnatives.pdf"), height = 4)
#jpeg(filename = file.path(figures, "Map_alldata.jpg"), quality = 100, res = 300, width = 2000, height = 2000)
mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF3, col="black", bg="black", cex= 1, pch=19)
points(dsSPDF4, col="red", bg="red", cex= 1, pch=19)
dev.off()

