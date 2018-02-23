## Libraries and functions
library(sp)
library(raster)
library(ape)
library(distrEx)
library(rgdal)
library(ncf)
library(Hmisc)
library(gplots)
library(distr)

#########################################################################################
## Variables in this section will need to be changed
#########################################################################################

source("C:\\Users\\hp39wasi\\sWorm\\SpatialGaps\\HellingersDistance\\HellsDistFunction.R")

#### Load Earthworm dataset
sites <- read.csv("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\3_Data\\Sites_2017-11-30.csv")
lat <- sites$Latitude__decimal_degrees
long <- sites$Longitude__Decimal_Degrees

#### Load Directory for Global Layers
GLdir <- "I:\\sDiv-postdocs-work\\Phillips\\sWorm\\SpatialGaps\\GlobalLayers\\rasters2\\"


#### Load global layers, for site level values only
precip <- raster(file.path(GLdir, "precip.tif"))
sites$precip <- extract(precip, data.frame(lat, long))
summary(sites$precip)

temp <- raster(file.path(GLdir, "temp.tif"))
sites$temp <- extract(temp, data.frame(lat, long))
summary(sites$temp)

p <- sites$precip
t <- sites$temp

#####################################################################################
## End of section - the next bit should run with no changes
#####################################################################################

###
## Method used in Gonzales et al 2016
quants <- c(1:4)
UD<- DiscreteDistribution(quants, c(0.25,0.25,0.25,0.25))

precip_HD <- list()
precip_HD[[1]] <- HellsDist(site_data_col = p, "precip.tif", inDir=GLdir, land.use = FALSE, n.reps = 1000, UD = UD)
precip_HD[[2]] <- HellsDistData(site_data_col = p,  "precip.tif", inDir = GLdir, land.use = FALSE, UD = UD)

temp_HD <- list()
temp_HD[[1]] <- HellsDist(site_data_col = t, "temp.tif", inDir=GLdir, land.use = FALSE, n.reps = 1000, UD = UD)
temp_HD[[2]] <- HellsDistData(site_data_col = t,  "temp.tif", inDir = GLdir, land.use = FALSE, UD = UD)


precip_HD_meansSD <- c(mean(precip_HD[[1]]), sd(precip_HD[[1]]))
precipdat_HD_meansSD <- c(precip_HD[[2]])

temp_HD_meansSD <- c(mean(temp_HD[[1]]), sd(temp_HD[[1]]))
tempdat_HD_meansSD <- c(temp_HD[[2]])

########### Plot
## pdf("C:\\Google Drive\\sWorm\\SpatialGaps-sWorm\\Figures\\HellingersDistance.pdf")
par(mar=c(3, 10, 1, 1))
pt.cex <- 1
pt.pch <- 19
layer_labels <- c("Precipitation Change", "Temperature Change")
plotCI(x=1, y=0, type="n", xlim=c(0,0.61), ylim=c(1,2), xlab = "Hellinger's d", yaxt = "n", ylab = "")
rect(-0.5,2.5, 1, 8,col="#A0A0A033",border=NA)

plotCI(x =  precip_HD_meansSD[1], y =2, uiw = precip_HD_meansSD[2],err='x', cex = pt.cex, gap =0, pch = pt.pch, add =TRUE)
plotCI(x =  precipdat_HD_meansSD[1], y =2,err='x',cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)

plotCI(x =  temp_HD_meansSD[1], y =1, uiw = temp_HD_meansSD[2],err='x', cex = pt.cex, gap =0, pch = pt.pch, add =TRUE)
plotCI(x =  tempdat_HD_meansSD[1], y =1,err='x',cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)
axis(2, at = 2:1, labels = layer_labels, las =2)
dev.off()



