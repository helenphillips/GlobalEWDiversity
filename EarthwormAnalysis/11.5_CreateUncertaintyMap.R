if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "TSGIS02"){
  setwd("C:/sWorm/EarthwormAnalysis")
}


if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\Users\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
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

source("Functions/cornerlabel2.R")

#######################################
# variables
######################################

wide_cm <- 12
wide_inch <- 4.75
point_size <- 10
plotlabcex <- 1
legendcex <- 0.8
resdpi <- 300
############################################################
## 'Global'
############################################################

bkg <- raster("I:\\sWorm\\ProcessedGLs-Dec2019\\Abundance\\CHELSA_bio10_1_AbundanceCutScaled.tif")
results_folder <- "C:\\Users\\hp39wasi\\Google Drive\\sWorm\\sWorm_WS\\EarthwormAnalysis\\Manuscripts\\SpatialAnalysis\\MapsFromZurich\\Earthworm_PCA_Chull_correction"

colbrks <-  c(0, seq(0.5, 1, length.out = 199))

r.cols <- viridis(199)


scale <- c(rep(r.cols[1], times = 20), rep(r.cols, each = 2), rep(r.cols[199], times = 20))

############################################################
## Species Richness
############################################################

richness1 <-  raster(file.path(results_folder, "Richness-0000000000-0000000000.tif"))
richness2 <-  raster(file.path(results_folder, "Richness-0000000000-0000032768.tif"))




png(file.path(figures, "Richness_CHull_Zurich_correction.png"),width=17.5,height=8.75,units="cm",res=resdpi)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))

par(mar=c(0.1,0.1,0.1,0.1))

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")

image(richness1, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(richness2, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")

corner.label2(label = "A", x = -1, y = 1, cex=plotlabcex, font = 2)

## Legend
par(mar=c(1,10,1,10))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )

mtext("50%", at =b[1], cex = legendcex)
mtext("100%", at = b[438], cex = legendcex)
mtext("Percentage of bands within interpolation range", side = 1, at = 250, cex = legendcex - 0.2)

dev.off()



############################################################
## BIOMASS
############################################################

biomass1 <-  raster(file.path(results_folder, "Biomass-0000000000-0000000000.tif"))
biomass2 <-  raster(file.path(results_folder, "Biomass-0000000000-0000032768.tif"))


png(file.path(figures, "Biomass_CHull_Zurich_correction.png"),width=17.5,height=8.75,units="cm",res=resdpi)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))

par(mar=c(0.1,0.1,0.1,0.1))

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")

image(biomass1, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(biomass2, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")

corner.label2(label = "C", x = -1, y = 1, cex=plotlabcex, font = 2)

## Legend
par(mar=c(1,10,1,10))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )


mtext("50%", at =b[1], cex = legendcex)
mtext("100%", at = b[438], cex = legendcex)
mtext("Percentage of bands within interpolation range", side = 1, at = 250, cex = legendcex - 0.2)

dev.off()




############################################################
## ABUNDANCE
############################################################


abundance1 <-  raster(file.path(results_folder, "Abundance-0000000000-0000000000.tif"))
abundance2 <-  raster(file.path(results_folder, "Abundance-0000000000-0000032768.tif"))


png(file.path(figures, "Abundance_CHull_Zurich_correction.png"),width=17.5,height=8.75,units="cm",res=resdpi)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))

par(mar=c(0.1,0.1,0.1,0.1))

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")

image(abundance1, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(abundance2, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")

corner.label2(label = "B", x = -1, y = 1, cex=plotlabcex, font = 2)

## Legend
par(mar=c(1,10,1,10))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )


mtext("50%", at =b[1], cex = legendcex)
mtext("100%", at = b[438], cex = legendcex)
mtext("Percentage of bands within interpolation range", side = 1, at = 250, cex = legendcex - 0.2)

dev.off()



