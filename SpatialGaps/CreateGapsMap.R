library(raster)
library(dismo)
library(maps)
setwd("C:/Users/hp39wasi/sWorm/SpatialGaps")

all_data <-"C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\3_Data"

files <- list.files(file.path(all_data))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
all_data <- read.csv(file.path(all_data, paste("sites_",date,".csv", sep = "")))

coord<-aggregate(cbind(all_data$Longitude__Decimal_Degrees, all_data$Latitude__decimal_degrees), list(all_data$Study_Name), mean)

coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")


# map <- raster("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\6_Data\\spRFinalRaster.tif")
mask <- raster("AccClimateMask.tif")

############################################
## Change in climate versus sampled points

#pdf(file = "ClimateChange_gaps.pdf")
jpeg(file = "ClimateChange_gaps.jpg", width = 1500, height = 1500)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
plot(mask, add = TRUE, legend = FALSE)
points(dsSPDF, col="black", bg="black", cex= 1.9, pch=19)
dev.off()

