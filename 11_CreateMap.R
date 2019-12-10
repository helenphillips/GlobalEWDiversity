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

# bkg <- raster("I:\\sWorm\\ProcessedGLs_revised\\CHELSA_bio10_1_BiomassCutScaled.tif")

#bkg <- raster("D:\\sWorm\\ProcessedGL\\Datasets\\bio10_1.tif")
bkg <- raster("I:\\sWorm\\ProcessedGLs-Dec2019\\Abundance\\CHELSA_bio10_1_AbundanceCutScaled.tif")

## The percentage that map is cut off at top and bottom

regions <- c("r12", "r13", "r14", "r15", "r16", "r17", "r18", "r19", 
             "r2", "r20", "r21", "r22", "r23", "r24", "r25", "r26", "r27", "r28",
             "r3", "r30", "r31", "r32", "r33", "r34", "r35", "r36", "r37", "r38", "r39",
             "r40", "r41", "r42", "r43", "r44", "r45", "r46", "r47", "r48", "r49",
             "r5", "r6", "r8", "r9")

############################################################
## Species Richness
############################################################

# results <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\Results\\Revised\\Richness"
# results <- "D:\\sWorm\\Results\\Richness"
results <- "I:\\sWorm\\Results\\Results-Dec2019\\Richness"


# regions <- c("africa", "asia", "europe", "latin_america", "north_america", "west_asia")

resultRaster <- "spRFinalRaster.tif"


minValues <- c()
maxValues <- c()

allValues <- c()

for(r in 1:length(regions)){
    
  print(r)
  tempR <- read.csv(file.path(results, regions[r], "predictedValues.csv"))
  
  if( !(all(is.na(tempR)))){
    
    minValues[r] <- min(tempR, na.rm = TRUE)
    maxValues[r] <- max(tempR, na.rm = TRUE)
    
    tempR <- na.omit(tempR)
    
    allValues <- c(allValues, as.vector(tempR[,1]))
  }else{print("All values are NAs")}
  

}

hist(allValues)
hist(allValues, xlim = c(0, 20), breaks = seq(0,  round(max(allValues, na.rm = TRUE)), by = 1))

# 
# hist(africa)
# hist(asia)
# hist(europe)
# hist(latin_america)
# hist(north_america)
# hist(west_asia)

minV <- min(minValues, na.rm = TRUE)
maxV <- max(maxValues, na.rm = TRUE)





colbrks <-  c(minV, seq(1, 4, length.out = 198), maxV)
# Changing the breaks
r.cols <- magma(199)

# png(file.path(figures, "Richness.png"),width=17.5,height=8.75,units="cm",res=resdpi)
pdf(file.path(figures, "Richness.pdf"),width= wide_inch, height= wide_inch/2, pointsize = point_size)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")


for(r in 1:length(regions)){
  tempR <- raster(file.path(results, regions[r], resultRaster))
  image(tempR, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
  
}

corner.label2(label = "B", x = -1, y = 1, cex=plotlabcex, font = 2)

## Legend
par(mar=c(0.1,13,1,13))
scale <- c(rep(magma(199)[1], times = 20), rep(magma(199), each = 2), rep(magma(199)[199], times = 20))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )
# b
# abline(v = 20, col = "white")
# abline(v = b[418], col = "black")

mtext(1, at = b[20], cex = legendcex)
mtext(6, at = b[418],cex = legendcex)
mtext("Number of species", at = 250, cex = legendcex)
dev.off()



############################################################
## BIOMASS
############################################################


# regions <- c("asia", "africa", "europe", "latin_america", "north_america", "west_asia")

# results <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\Results\\Revised\\Biomass"
# results <- "D:\\sWorm\\Results\\Biomass"

results <- "I:\\sWorm\\Results\\Results-Dec2019\\Biomass"


resultRaster <- "BiomassFinalRaster.tif"


minValues <- c()
maxValues <- c()

allValues <- c()

for(r in 1:length(regions)){
  
  print(r)
  tempR <- read.csv(file.path(results, regions[r], "predictedValues.csv"))
  
  #back transform!
  minValues[r] <- min(tempR, na.rm = TRUE)
  maxValues[r] <- max(tempR, na.rm = TRUE)
  
  tempR <- na.omit(tempR)
  
  allValues <- c(allValues, as.vector(tempR))
  
}

hist(allValues)
hist(allValues, xlim = c(0, 20), breaks = seq(0,  round(max(allValues, na.rm = TRUE)), by = 1))


#   
# hist(africa, ylim =c(1, 1000))
# hist(asia, ylim =c(1, 1000))
# hist(europe, xlim = c(0, 50), breaks = seq(minValue(europe),  maxValue(europe), by = 1))
# hist(latin_america, ylim =c(1, 1000))
# hist(north_america, ylim =c(1, 1000))
# hist(west_asia, ylim =c(1, 1000))
 minV <- min(minValues, na.rm = TRUE)
 maxV <- max(maxValues, na.rm = TRUE)
 
#colbrks <-  c(minV, seq(1, 150, length.out = 198), maxV)
colbrks <-  c(minV, seq(-3, 15, length.out = 198), maxV)

r.cols <- magma(199)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))


#png(file.path(figures, "Biomass.png"),width=17.5,height=8.75,units="cm",res=resdpi)
pdf(file.path(figures, "Biomass.pdf"),width= wide_inch, height= wide_inch/2, pointsize = point_size)

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")

for(r in 1:length(regions)){
  tempR <- raster(file.path(results, regions[r], resultRaster))
  image(tempR, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
  
}

corner.label2(label = "D", x = -1, y = 1, cex = plotlabcex, font = 2)


## Legend
par(mar=c(1,13,1,13))
scale <- c(rep(magma(199)[1], times = 20), rep(magma(199), each = 2), rep(magma(199)[199], times = 20))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )
mtext("1", at =b[20], cex = legendcex)
mtext("150", at = b[418], cex = legendcex)
mtext(expression("Biomass (grams per" ~ m^{2} ~ ")"), at = 250, cex = legendcex)
dev.off()



############################################################
## ABUNDANCE
############################################################




# results <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\Results\\Revised\\Abundance"
# results <- "D:\\sWorm\\Results\\Abundance"
results <- "I:\\sWorm\\Results\\Results-Dec2019\\Abundance"

resultRaster <- "AbundanceFinalRaster.tif"


minValues <- c()
maxValues <- c()

allValues <- c()


for(r in 1:length(regions)){
  
  print(r)
  tempR <- read.csv(file.path(results, regions[r], "predictedValues.csv"))
  
  #back transform!
  minValues[r] <- min(tempR, na.rm = TRUE)
  maxValues[r] <- max(tempR, na.rm = TRUE)
  
  tempR <- na.omit(tempR)
  
  allValues <- c(allValues, as.vector(tempR))
  
}

allValues <- exp(allValues ) - 1
hist(allValues)
hist(allValues, xlim = c(0, 500), breaks = seq(floor(min(allValues, na.rm = TRUE)),  ceiling(max(allValues, na.rm = TRUE)), by = 1))


minV <- min(minValues, na.rm = TRUE)
maxV <- max(maxValues, na.rm = TRUE)
# hist(asia)
# hist(africa)
# hist(europe)
# hist(north_america)
# hist(latin_america)
# hist(west_asia)



###########################################

colbrks <-  c(minV, seq(5, 100, length.out = 198), maxV)

r.cols <- magma(199)

#png(file.path(figures, "Abundance.png"),width=17.5,height=8.75,units="cm",res=resdpi)
pdf(file.path(figures, "Abundance.pdf"),width= wide_inch, height= wide_inch/2, pointsize = point_size)

nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")

for(r in 1:length(regions)){
  tempR <- raster(file.path(results, regions[r], resultRaster))
  tempR <- exp(tempR) - 1
  image(tempR, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
  
}
corner.label2(label = "C", x = -1, y = 1, cex = plotlabcex, font = 2)


## Legend
par(mar=c(1,13,1,13))
scale <- c(rep(magma(199)[1], times = 20), rep(magma(199), each = 2), rep(magma(199)[199], times = 20))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )
mtext("5", at = b[20], cex = legendcex)
mtext("150", at = b[418], cex = legendcex)
mtext(expression("Abundance (individuals per" ~ m^{2} ~ ")"), at = 250, cex = legendcex)

dev.off()



