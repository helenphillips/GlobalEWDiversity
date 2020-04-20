

## Create map
library(raster)
library(RColorBrewer)
library(maps)
library(mapdata)
library(maptools)
library(sp)
library(viridis)


if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "TSGIS02"){
  setwd("C:/sWorm/EarthwormAnalysis")
  rasterOptions(tmpdir = "D:\\sWorm\\trash", chunksize = 524288, maxmemory = 134217728)
  
}


if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\Users\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
  
}


figures <- "Figures"
source("Functions/cornerlabel2.R")

#######################################
# variables
######################################

newRasters <- "D:\\sWorm\\scaledRasters\\new"
oldRasters <- "D:\\sWorm\\scaledRasters\\old"

############################################################
## 'Global'
############################################################

## The percentage that map is cut off at top and bottom


############################################################
## New Species Richness
############################################################

regions <- c("r12", "r13", "r14", "r15", "r16", "r17", "r18", "r19", 
             "r2", "r20", "r21", "r22", "r23", "r24", "r25", "r26", "r27", "r28",
             "r3", "r30", "r31", "r32", "r33", "r34", "r35", "r36", "r37", "r38", "r39",
             "r40", "r41", "r42", "r43", "r44", "r45", "r46", "r47", "r48", "r49",
             "r5", "r6", "r8", "r9")

results <- "I:\\sWorm\\Results\\Results-Dec2019\\Richness"


resultRaster <- "spRFinalRaster.tif"


minValues <- c()
maxValues <- c()

allValues <- c()

for(r in 1:length(regions)){
  
  print(r)
  tempR <- read.csv(file.path(results, regions[r], "predictedValues.csv"), header = FALSE)
  

  
  if( !(all(is.na(tempR)))){
    
    minValues[r] <- min(tempR, na.rm = TRUE)
    maxValues[r] <- max(tempR, na.rm = TRUE)
    
    tempR <- na.omit(tempR)
    
    allValues <- c(allValues, as.vector(tempR[,1]))
  }else{print("All values are NAs")}
  
  
}


minV <- min(minValues, na.rm = TRUE)
maxV <- max(maxValues, na.rm = TRUE)


# allValues <- (allValues - minV) / (maxV - minV)


for(r in 1:length(regions)){
  tempR <- raster(file.path(results, regions[r], resultRaster))
  

  
  tempR <- (tempR - minV) / (maxV - minV)
  tempR <- writeRaster(tempR, filename=file.path(newRasters, paste0(regions[r], ".tif")), format="GTiff", overwrite=TRUE)
  
  
}


############################################################
## old Species Richness
############################################################

regions <- c("asia", "africa", "europe", "west_asia", "north_america", "latin_america")

results <- "I:\\sWorm\\Results\\Richness"


resultRaster <- "spRFinalRaster.tif"


minValues <- c()
maxValues <- c()

allValues <- c()

for(r in 1:length(regions)){
  
  print(r)
  tempR <- read.csv(file.path(results, regions[r], "predictedValues.csv"), header = FALSE)
  
  tempR <- exp(tempR)
  
  if( !(all(is.na(tempR)))){
    
    minValues[r] <- min(tempR, na.rm = TRUE)
    maxValues[r] <- max(tempR, na.rm = TRUE)
    
    tempR <- na.omit(tempR)
    
    allValues <- c(allValues, as.vector(tempR[,1]))
  }else{print("All values are NAs")}
  
  
}


minV <- min(minValues, na.rm = TRUE)
maxV <- max(maxValues, na.rm = TRUE)


# allValues <- (allValues - minV) / (maxV - minV)


for(r in 1:length(regions)){
  tempR <- raster(file.path(results, regions[r], resultRaster))
  tempR <- exp(tempR)
  
  tempR <- (tempR - minV) / (maxV - minV)
  tempR <- writeRaster(tempR, filename=file.path(oldRasters, paste0(regions[r], ".tif")), format="GTiff", overwrite=TRUE)
  
  
}


##########################
## PLOTTING
############################

diff <- raster("D:\\sWorm\\scaledRasters\\diff")

# minV <- cellStats(diff, min, na.rm = TRUE)
# maxV <- cellStats(diff, max, na.rm = TRUE)

minV <- -0.8980593
maxV <- 0.9781621

colbrks <-  c(minV, -0.001, 0.001,maxV)

# Changing the breaks
r.cols <- magma(3)


png(file.path(figures, "Richness_scaledDiff_001.png"),width=17.5,height=8.75,units="cm",res=resdpi)

nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))
image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
image(diff, ylim = c(-90, 90), xlim = c(-180, 180), add = TRUE, col=r.cols, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="", bty="n")


corner.label2(label = "0.001", x = -1, y = 1, cex=plotlabcex, font = 2)

## Legend
par(mar=c(1,10,1,10))
scale <- c(rep(magma(3)[1], times = 146), rep(magma(3)[2], each = 146), rep(magma(3)[3], times = 146))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )


mtext("Old is lower", at = b[100], cex = legendcex)
mtext("Same", at = b[210], cex = legendcex)
mtext("Old is higher", at = b[350],cex = legendcex)


dev.off()

colbrks <-  c(minV, -0.01, 0.01,maxV)

# Changing the breaks
r.cols <- magma(3)


png(file.path(figures, "Richness_scaledDiff_01.png"),width=17.5,height=8.75,units="cm",res=resdpi)

nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))
image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
image(diff, ylim = c(-90, 90), xlim = c(-180, 180), add = TRUE, col=r.cols, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="", bty="n")


corner.label2(label = "0.01", x = -1, y = 1, cex=plotlabcex, font = 2)

## Legend
par(mar=c(1,10,1,10))
scale <- c(rep(magma(3)[1], times = 146), rep(magma(3)[2], each = 146), rep(magma(3)[3], times = 146))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )


mtext("Old is lower", at = b[100], cex = legendcex)
mtext("Same", at = b[210], cex = legendcex)
mtext("Old is higher", at = b[350],cex = legendcex)


dev.off()
