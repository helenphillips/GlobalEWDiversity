if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}


if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

## Create map
library(raster)
library(RColorBrewer)
library(maps)
library(mapdata)
library(maptools)
library(sp)
library(viridis)

############################################################
## 'Global'
############################################################

bkg <- raster("I:\\sWorm\\ProcessedGLs\\CHELSA_bio10_1_BiomassCutScaled.tif")

## The percentage that map is cut off at top and bottom




############################################################
## Species Richness
############################################################

results <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\Results\\Richness"

# regions <- c("africa", "asia", "europe", "latin_america", "north_america", "west_asia")

resultRaster <- "spRFinalRaster.tif"

africa <- raster(file.path(results, "africa", resultRaster))
africa <- exp(africa)
asia <-  raster(file.path(results, "asia", resultRaster))
asia <- exp(asia)
europe <- raster(file.path(results, "europe", resultRaster))
europe <- exp(europe) 
latin_america <- raster(file.path(results, "latin_america", resultRaster))
latin_america <- exp(latin_america) 
north_america <- raster(file.path(results, "north_america", resultRaster))
north_america <- exp(north_america) 
west_asia <- raster(file.path(results, "west_asia", resultRaster))
west_asia <- exp(west_asia)

hist(africa)
hist(asia)
hist(europe)
hist(latin_america)
hist(north_america)
hist(west_asia)


minV <-min(c(minValue(africa),
             minValue(asia),
             minValue(europe),
             minValue(latin_america),
             minValue(north_america),
             minValue(west_asia)))

maxV <-max(c(maxValue(africa),
             maxValue(asia),
             maxValue(europe),
             maxValue(latin_america),
             maxValue(north_america),
             maxValue(west_asia)))

meanV <- mean(c(cellStats(africa, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(asia, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(europe, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(latin_america, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(north_america, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(west_asia, stat='mean', na.rm=TRUE, asSample=TRUE)))

sdV <- mean(c(cellStats(africa, stat='sd', na.rm=TRUE, asSample=TRUE),
                cellStats(asia, stat='sd', na.rm=TRUE, asSample=TRUE),
                cellStats(europe, stat='sd', na.rm=TRUE, asSample=TRUE),
                cellStats(latin_america, stat='sd', na.rm=TRUE, asSample=TRUE),
                cellStats(north_america, stat='sd', na.rm=TRUE, asSample=TRUE),
                cellStats(west_asia, stat='sd', na.rm=TRUE, asSample=TRUE)))

# breakpoint_top <- 0.3
# breakpoint_bottom <- 0.02
# diff <- maxV - minV
# 
# top20 <- maxV - (diff * breakpoint_top)
# bottom20 <- minV + (diff * breakpoint_bottom)


colbrks <-  c(minV, seq(1, 6, length.out = 198), maxV)

# seq(minV, maxV, length.out = 200)
# actual <- c('#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
#            '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')
# r.cols <- colorRampPalette(actual, space = "rgb")




r.cols <- magma(199)

png(file.path(figures, "Richness.png"),width=17.5,height=8.75,units="cm",res=1200)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="")


image(africa, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(asia, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(europe, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(latin_america, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(north_america, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(west_asia, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")

## Legend
par(mar=c(1,13,1,13))
scale <- c(rep(magma(199)[1], times = 20), rep(magma(199), each = 2), rep(magma(199)[199], times = 20))
barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )
mtext(1, at = 75, cex = 1)
mtext(6, at = 460, cex = 1)
mtext("Number of species", at = 250, cex = 0.5)
dev.off()



############################################################
## BIOMASS
############################################################


# regions <- c("africa", "asia", "europe", "latin_america", "north_america", "west_asia")

results <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\Results\\Biomass"

 regions <- c("asia", "africa", "europe", "latin_america", "north_america", "west_asia")
 resultRaster <- "BiomassFinalRaster.tif"
 
 africa <- raster(file.path(results, "africa", resultRaster))
 # africa <- exp(africa) - 1
 asia <-  raster(file.path(results, "asia", resultRaster))
 # asia <- exp(asia) - 1
 europe <- raster(file.path(results, "europe", resultRaster))
 # europe <- exp(europe) - 1
 latin_america <- raster(file.path(results, "latin_america", resultRaster))
 # latin_america <- exp(latin_america) - 1
 north_america <- raster(file.path(results, "north_america", resultRaster))
 # north_america <- exp(north_america) - 1
 west_asia <- raster(file.path(results, "west_asia", resultRaster))
 # west_asia <- exp(west_asia) - 1
 
 
 hist(africa, ylim =c(1, 1000))
 hist(asia, ylim =c(1, 1000))
 hist(europe, ylim =c(1, 1000))
 hist(latin_america, ylim =c(1, 1000))
 hist(north_america, ylim =c(1, 1000))
 hist(west_asia, ylim =c(1, 1000))
 
minV <-min(c(minValue(africa),
             minValue(asia),
             minValue(europe),
             minValue(latin_america),
             minValue(north_america),
             minValue(west_asia)))
 
maxV <-max(c(maxValue(africa),
             maxValue(asia),
             maxValue(europe),
             maxValue(latin_america),
             maxValue(north_america),
             maxValue(west_asia)))

meanV <- mean(c(cellStats(africa, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(asia, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(europe, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(latin_america, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(north_america, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(west_asia, stat='mean', na.rm=TRUE, asSample=TRUE)))

sdV <- mean(c(cellStats(africa, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(asia, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(europe, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(latin_america, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(north_america, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(west_asia, stat='sd', na.rm=TRUE, asSample=TRUE)))
# 
 

colbrks <-  c(minV, seq(1, 6, length.out = 198), maxV)

r.cols <- magma(199)

png(file.path(figures, "Biomass.png"),width=17.5,height=8.75,units="cm",res=1200)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))
image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="")

image(africa, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(asia, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(europe, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(latin_america, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(north_america, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(west_asia, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")

## Legend
par(mar=c(1,13,1,13))
scale <- c(rep(magma(199)[1], times = 20), rep(magma(199), each = 2), rep(magma(199)[199], times = 20))
barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )
mtext("1", at = 75, cex = 1)
mtext("6", at = 460, cex = 1)
mtext("log-Biomass (g/m2)", at = 250, cex = 0.5)
dev.off()



############################################################
## ABUNDANCE
############################################################




results <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\Results\\Abundance"

regions <- c("asia", "afria", "europe", "latin_america", "north_america", "west_asia")
resultRaster <- "AbundanceFinalRaster.tif"

africa <- raster(file.path(results, "africa", resultRaster))
africa <- exp(africa) - 1
asia <-  raster(file.path(results, "asia", resultRaster))
asia <- exp(asia) - 1
europe <- raster(file.path(results, "europe", resultRaster))
europe <- exp(europe) - 1
latin_america <- raster(file.path(results, "latin_america", resultRaster))
latin_america <- exp(latin_america) - 1
north_america <- raster(file.path(results, "north_america", resultRaster))
north_america <- exp(north_america) - 1
west_asia <- raster(file.path(results, "west_asia", resultRaster))
west_asia <- exp(west_asia) - 1

minV <-min(c(minValue(africa),
             minValue(asia),
             minValue(europe),
             minValue(latin_america),
             minValue(north_america),
             minValue(west_asia)))

maxV <-max(c(maxValue(africa),
             maxValue(asia),
             maxValue(europe),
             maxValue(latin_america),
             maxValue(north_america),
             maxValue(west_asia)))

hist(asia)
hist(africa)
hist(europe)
hist(north_america)
hist(latin_america)
hist(west_asia)

meanV <- mean(c(cellStats(africa, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(asia, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(europe, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(latin_america, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(north_america, stat='mean', na.rm=TRUE, asSample=TRUE),
                cellStats(west_asia, stat='mean', na.rm=TRUE, asSample=TRUE)))

sdV <- mean(c(cellStats(africa, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(asia, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(europe, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(latin_america, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(north_america, stat='sd', na.rm=TRUE, asSample=TRUE),
              cellStats(west_asia, stat='sd', na.rm=TRUE, asSample=TRUE)))



###########################################

colbrks <-  c(minV, seq(5, 150, length.out = 198), maxV)

r.cols <- magma(199)

png(file.path(figures, "Abundance.png"),width=17.5,height=8.75,units="cm",res=1200)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="")

image(africa, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(asia, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(europe, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(latin_america, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(north_america, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(west_asia, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")

## Legend
par(mar=c(1,13,1,13))
scale <- c(rep(magma(199)[1], times = 20), rep(magma(199), each = 2), rep(magma(199)[199], times = 20))
barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )
mtext("5 ind/m2", at = 75, cex = 1)
mtext("150 ind/m2", at = 460, cex = 1)
mtext("Abundance", at = 250, cex = 0.5)

dev.off()



