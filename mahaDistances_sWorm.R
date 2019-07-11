## Computation of the distances

# The mahnalobis distance gives an indication of how far is a given point from the distribution.
# Create a data frame that only contains the predictor variables.
# This data frame is called predictorsData

# setwd("C:/Research/Projects/sWorm/Files")

library(raster)
library(sp)
library(rgdal)
# Load you data here:
# data = read.csv("richness.csv",header=TRUE)
data <- read.csv("temp\\sitesRichness_2019-07-03.csv")
summary(data)

# Removal of no data (this is only valuable if you populate your dataset directly from the maps)
# data[data == -9999] <- NA
#data = na.omit(data)
# summary(data)

# Create a dataset in a given order - Land Cover is not included here because we are working with averages and LC is a categorical variable
predictorsData = data.frame(ID = data$ID,
                       #     bio10_1 = (data$CHELSA_bio10_1),
#bio10_4  = (data$CHELSA_bio10_4),
bio10_7  = (data$CHELSA_bio10_7),
#bio10_12  = (data$CHELSA_bio10_12),
bio10_15  = (data$CHELSA_bio10_15),
cation = (data$CECSOL_weighted),
elevation  = (data$elevation),
ai           = (data$ai_yr_TIF),
ph           = (data$PHIHOX_weighted),
pet          = (data$pet_he_yr_TIF),
#pet_sd       = (data$pet_he_SD),
#snow         = (data$snow_2015_sum),
clay         = (data$CLYPPT_weighted),
silt         = (data$SLTPPT_weighted),
carbon       = (data$ORCDRC_weighted))

summary(predictorsData[,-1])


# These two lines compute the covariance matrix and the vector of means
# of the predictors. This are a 15x15 matrix and a vector with 15 entries.

# TODO : THIS SHOULD BE DONE AFTER THE DATA IS SCALED

sigma<-cov(predictorsData[,-1])
mu<-colMeans(predictorsData[,-1]) # in both cases take care what columns are you eliminating (IDs and Coordinates should not be here)


## Now need to cap the global layers
# And create a dataframe of all the values


scaleGL <- function(layername, data, load = ".", save="."){
  if(!(layername %in% names(data))) stop("Layername must be in data")
  print(layername)
  
  dat <- data[,which(names(data) == layername)]
  
  min <- min(dat, na.rm = TRUE) #NOT SCALED YET
  max <- max(dat, na.rm = TRUE)
  
  tif <- raster(load)

  print("Capping minimum value")
  tif[tif<min] <- min
  print("Capping maximum value")
  tif[tif>max] <- max
  
  print(cellStats(tif, stat='mean', na.rm=TRUE, asSample=TRUE))
  print(cellStats(tif, stat='sd', na.rm=TRUE, asSample=TRUE))
  
  
  dat_mean <- mean(dat)
  print(dat_mean)
  dat_sd <- sd(dat)
  print(dat_sd)
  
  
  # print("Scaling tif")
  # tif <- (tif - dat_mean)/dat_sd
  
  # print("Saving tif")
  writeRaster(tif, save, format = "GTiff")
  # do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))
  
   rm(tif)
}

savefolder <- "D:\\sWorm\\ProcessedGL\\10km_cutscaled"
folder <- "D:\\sWorm\\ProcessedGL\\Datasets\\Same_resolution_10k"
scaleGL(layername = "carbon", data = predictorsData, load =  file.path(folder, "carbon.tif"),
        save = file.path(savefolder, "carbon.tif"))
scaleGL(layername = "ph", data = predictorsData, load =  file.path(folder, "ph.tif"),
        save = file.path(savefolder, 'ph.tif'))
scaleGL(layername = "clay", data = predictorsData, load =  file.path(folder, "clay.tif"),
        save = file.path(savefolder, 'clay.tif'))
scaleGL(layername = "silt", data = predictorsData, load =  file.path(folder, "silt.tif"),
        save = file.path(savefolder, 'silt.tif'))
scaleGL(layername = "cation", data = predictorsData, load =  file.path(folder, "cation.tif"),
        save = file.path(savefolder, 'cation.tif'))

scaleGL(layername = 'bio10_1', data = predictorsData, load = file.path(folder, 'bio10_1.tif'),
        save = file.path(savefolder, 'bio10_1.tif'))
scaleGL(layername = 'bio10_4', data = predictorsData, load = file.path(folder, 'bio10_4.tif'),
        save = file.path(savefolder, 'bio10_4.tif'))
scaleGL(layername = 'bio10_7', data = predictorsData, load = file.path(folder, 'bio10_7.tif'),
        save = file.path(savefolder, 'bio10_7.tif'))

scaleGL(layername = 'bio10_12', data = predictorsData, load = file.path(folder, 'bio10_12.tif'),
        save = file.path(savefolder, 'bio10_12.tif'))
scaleGL(layername = 'bio10_15', data = predictorsData, load = file.path(folder, 'bio10_15.tif'),
        save = file.path(savefolder, 'bio10_15.tif'))

scaleGL(layername = 'ai', data = predictorsData, load = file.path(folder, 'ai.tif'),
        save = file.path(savefolder, 'ai.tif'))
scaleGL(layername = 'pet', data = predictorsData, load = file.path(folder, 'pet.tif'),
        save = file.path(savefolder, 'pet.tif'))
scaleGL(layername = 'pet_sd', data = predictorsData, load = file.path(folder, 'pet_sd.tif'),
        save = file.path(savefolder, 'pet_sd.tif'))
scaleGL(layername = 'elevation', data = predictorsData, load = file.path(folder, 'elevation.tif'),
        save = file.path(savefolder, 'elevation.tif'))

# Read the file with the new input variables - To be mapped - The shapefile that I am going to send you
distance=readOGR("I:\\sWorm\\Dataset_072019", layer = "point_grid") 
# distance <- read.csv("Dataset_072019\\data_globe_v2.txt",header=TRUE)

distance <- data.frame(distance@data)

summary(distance)



tifs <- c("ph", "carbon", "cation", "silt", "clay",  "bio10_7",  "bio10_15",
          "ai", "pet",  "elevation")


distanceData <- data.frame(pointID = distance$POINTID, X = distance$X, Y = distance$Y)
for(t in tifs){
  tif <- raster(file.path(savefolder, 
                          paste(t, ".tif", sep="")))
  distanceData$V1 <- extract(tif, data.frame(distanceData$X, distanceData$Y))
  names(distanceData)[names(distanceData) == "V1"] <- t
}


for(t in tifs){
  print(tif)
  tif <- raster(file.path(savefolder, 
                          paste(t, ".tif", sep="")))
 print(cellStats(tif, stat = "mean"))
 print(cellStats(tif, stat = "max"))
 print(cellStats(tif, stat = "min"))
}



distanceData <- distanceData[,c(1,2, 3, 9, 10, 6, 13, 11, 4, 12, 8, 7, 5)]

# 
# # Create a new dataframe - keep the same order
# distanceData = with(data, data.frame(distance$POINTID))
# distanceData$bio10_1      <- distance$bio10_1
# distanceData$bio10_4      <- distance$bio10_4
# distanceData$bio10_7      <- distance$bio10_7
# distanceData$bio10_12     <- distance$bio10_12
# distanceData$bio10_15     <- distance$bio10_15
# distanceData$cation       <- distance$cation
# distanceData$elevation    <- distance$elevation
# distanceData$ai           <- distance$ai
# distanceData$ph           <- distance$ph
# distanceData$pet          <- distance$pet
# distanceData$pet_sd       <- distance$pet_sd
# distanceData$snow         <- distance$snow
# distanceData$clay         <- distance$clay
# distanceData$silt         <- distance$silt
# distanceData$carbon       <- distance$carbon

# Remove NAs if existing  (this method does not cope with NAs)
# distanceData[distanceData == -9999] <- NA
distanceData = distanceData[complete.cases(distanceData),]
summary(distanceData[,-1])
d = (distanceData[,-c(1:3)]) # create a subset dataset to be easier to handle - not critical and you can remove it if you prefer

# Create a new column in the dataframe that saves the distance for each point.
distanceData$mahaDistance <- mahalanobis(d,mu,sigma) # Calculate the multidimentional distance
summary(distanceData)



length(which(distanceData$mahaDistance < qchisq(.975, df=10))) # How many points are in the 97.5% quantile
# 465034
100 * (length(which(distanceData$mahaDistance < qchisq(.975, df=10))) / length(distanceData$pointID)) # estimating the proportion of pixels included (please note that the statistical definition of an outlier here is .975)
# 32.19543
distanceData_final <- data.frame(pointID = distanceData$pointID, X= distanceData$X,Y= distanceData$Y)
distanceData_final$mahaDist = distanceData$mahaDistance

write.csv(distanceData_final,"mahaDistances_richness_v1.csv") # Write the file with the distance.

###############

r <- raster(file.path(folder, "ph.tif"))
dimensions <- dim(r)
resol <-res(r)
coordred <- crs(r)
exten <- extent(r)


# CSV to shapefile


coordinates(distanceData_final)=~X+Y
proj4string(distanceData_final)<- coordred
DD <-spTransform(distanceData_final,CRS("+proj=longlat"))

raster::shapefile(DD, "MahalaDistance.shp")


rast <- raster()
extent(rast) <- extent(DD) # this might be unnecessary
ncol(rast) <- dim(r)[2] # this is one way of assigning cell size / resolution
nrow(rast) <- dim(r)[1]

# And then ... rasterize it! This creates a grid version 
# of your points using the cells of rast, values from the IP field:
rast2 <- rasterize(DD, rast, DD$mahaDist, fun=mean) 

image(rast2)

## nic plot

library(RColorBrewer)
library(maps)
library(mapdata)
library(maptools)
library(sp)
library(viridis)


cutoff <- qchisq(.975, df=10)
minV <- min(distanceData$mahaDistance)
maxV <- max(distanceData$mahaDistance)
colbrks <-  c(seq( minV, maxV, length = 199))

r.cols <- viridis(198)
pdf(file.path("Figures", "Richness_MD.pdf"),width=(17.5/2),height=(8.75/2))
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))

image(rast2, col=r.cols, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")

# plotrix::corner.label(label = "(b)", x = -1, y = 1, cex=0.1)

## Legend
par(mar=c(1,13,1,13))
scale <- c(rep(viridis(198), each = 4))
b <- barplot(rep(1, (198 * 4)), col = scale, border =scale, axes = FALSE )
# b
# abline(v = 20, col = "white")
# abline(v = b[418], col = "black")

mtext(round(minV, digits = 2), at = b[1], cex = 1)
mtext(round(maxV, digits = 2), at = b[792], cex = 1)
mtext("M-Distance", at = 250, cex = 0.5)



cutoff

198 = 233 # (maxV - minV)
198 / 233 = 1
0.84 = 1


cutoff - minV 
14.66 
0.84 * 14.66

12 * 4

mtext(round(cutoff, digits = 2), at = b[48], cex = 1)



cutoff - minV # 14.66



dev.off()





##### 

europe <- distanceData[distanceData$Y > 30,]
europe <- europe[europe$X > -16,]
europe <- europe[europe$X < 50,]


