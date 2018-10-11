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

############################################################
## Functional Richness
############################################################

results <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\Results\\FunctionalRichness"
resultRaster <- "FGRichnessFinalRaster.tif"

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
colbrks <-  c(minV, seq(0.5,3, lengthou.out = 198), maxV)

r.cols <- magma(199)

png(file.path(figures, "FunctionalRichness.png"),width=17.5,height=8.75,units="cm",res=1200)
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
mtext("0.5", at = 75, cex = 1)
mtext("3", at = 430, cex = 1)
mtext("Functional Richness", at = 250, cex = 0.5)
dev.off()


#####################
## A discreet map

colbrks <-  c(minV, seq(0.5, 4.5, by = 1), maxV)

r.cols <- magma(6)

png(file.path(figures, "FunctionalRichness_discreet.png"),width=17.5,height=8.75,units="cm",res=1200)
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
scale <- rep(magma(6), each = 20)
barplot(rep(1, 120), col = scale, border =scale, axes = FALSE )
mtext("0", at = 10, cex = 1)
mtext("1", at = 35, cex = 1)
mtext("2", at = 60, cex = 1)
mtext("3", at = 85, cex = 1)
mtext("4", at = 110, cex = 1)
mtext(">4", at = 130, cex = 1)
mtext("Functional Richness", at = 80, cex = 0.5, line = -1.5)
dev.off()


