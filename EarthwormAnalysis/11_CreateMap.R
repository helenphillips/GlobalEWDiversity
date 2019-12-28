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
legendcex <- 0.9
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
  tempR <- read.csv(file.path(results, regions[r], "predictedValues.csv"), header = FALSE)
  
  if( !(all(is.na(tempR)))){
    
    minValues[r] <- min(tempR, na.rm = TRUE)
    maxValues[r] <- max(tempR, na.rm = TRUE)
    
    tempR <- na.omit(tempR)
    
    allValues <- c(allValues, as.vector(tempR[,1]))
  }else{print("All values are NAs")}
  

}
quantile(allValues, probs = seq(0.98, 1, by = 0.001))
# 99% of the data is below 4 (species)
# 98% of the data is below 3 species

hist(allValues)
hist(allValues, xlim = c(0, 20), breaks = seq(0,  round(max(allValues, na.rm = TRUE)), by = 1))


minV <- min(minValues, na.rm = TRUE)
maxV <- max(maxValues, na.rm = TRUE)

minV <- log(minV)
maxV <- log(maxV)

colbrks <-  c(minV, seq(log(1), log(4), length.out = 198), maxV)
# Changing the breaks
 r.cols <- magma(199)

# r.cols <- plasma(199)


# png(file.path(figures, "Richness.png"),width=17.5,height=8.75,units="cm",res=resdpi)
# pdf(file.path(figures, "Richness.pdf"),width= wide_inch, height= wide_inch/2, pointsize = point_size)
# pdf(file.path(figures, "Richness_noLabel.pdf"),width= wide_inch, height= wide_inch/2, pointsize = point_size)
#pdf(file.path(figures, "Richness_correction.pdf"),width= wide_inch, height= wide_inch/2, pointsize = point_size)
# png(file.path(figures, "Richness_correction.png"),width=17.5,height=8.75,units="cm",res=resdpi)
png(file.path(figures, "Richness_correction_logged.png"),width=17.5,height=8.75,units="cm",res=resdpi)

nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")


for(r in 1:length(regions)){
  tempR <- raster(file.path(results, regions[r], resultRaster))
  tempR <- log(tempR)
  image(tempR, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
  
}

# corner.label2(label = "B", x = -1, y = 1, cex=plotlabcex, font = 2)

## Legend
par(mar=c(1,10,1,10))
scale <- c(rep(magma(199)[1], times = 20), rep(magma(199), each = 2), rep(magma(199)[199], times = 20))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )
# b
# abline(v = 20, col = "white")
# abline(v = b[418], col = "black")

mtext(1, at = b[20], cex = legendcex)
mtext(4, at = b[418],cex = legendcex)
mtext("Number of species", side = 1, at = 250, cex = legendcex - 0.2)
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
  tempR <- read.csv(file.path(results, regions[r], "predictedValues.csv"), header = FALSE)
  
  
  if( !(all(is.na(tempR)))){
    tempR[,1] <- exp(tempR[,1]) - 1
    
    minValues[r] <- min(tempR, na.rm = TRUE)
    maxValues[r] <- max(tempR, na.rm = TRUE)
    
    tempR <- na.omit(tempR)
    
    allValues <- c(allValues, as.vector(tempR[,1]))
  }else{print("All values are NAs")}
  
  
}

mean(allValues)
median(allValues)
sd(allValues)

hist(allValues)
hist(allValues, xlim = c(0, 100),  breaks = seq( floor(min(minValues, na.rm = TRUE)),  round(max(allValues, na.rm = TRUE) + 1), by = 1))

quantile(allValues, probs = seq(0.96, 1, by = 0.001))
quantile(allValues, probs = seq(0.99, 1, by = 0.0001))

## For BIOMASS

unrealistic <-2000 # i.e. an unrealistic sample of earthworms is 2kg
## In Lavelle's book (soil Ecology), he said abundance could be 2000ind m-2
## Assuming 1g per individual 
## His biomass estimates are max 160g m-2

(length(which(allValues > unrealistic ))) / length(allValues) * 100
# 0.1232064
 

minV <- min(minValues, na.rm = TRUE)
 maxV <- max(maxValues, na.rm = TRUE)
 
colbrks <-  c(minV, seq(1, 150, length.out = 198), maxV)

r.cols <- magma(199)



#png(file.path(figures, "Biomass.png"),width=17.5,height=8.75,units="cm",res=resdpi)
#pdf(file.path(figures, "Biomass.pdf"),width= wide_inch, height= wide_inch/2, pointsize = point_size)
png(file.path(figures, "Biomass_correction.png"),width=17.5,height=8.75,units="cm",res=resdpi)

nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
# layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))

image(bkg, ylim = c(-90, 90), xlim = c(-180, 180), col = "gray90", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")

for(r in 1:length(regions)){
  tempR <- raster(file.path(results, regions[r], resultRaster))
  tempR <- exp(tempR) - 1
  image(tempR, col=r.cols, add = TRUE, breaks=colbrks, xaxt="n", yaxt="n", ylab="", xlab="")
  
}

corner.label2(label = "D", x = -1, y = 1, cex = plotlabcex, font = 2)


## Legend
par(mar=c(1,10,1,10))
scale <- c(rep(magma(199)[1], times = 20), rep(magma(199), each = 2), rep(magma(199)[199], times = 20))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )
mtext("1", at =b[20], cex = legendcex)
mtext("150", at = b[418], cex = legendcex)
mtext(expression("Biomass (grams per" ~ m^{2} ~ ")"), side = 1, at = 250, cex = legendcex - 0.2)
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
#  for(r in 1:2){
  print(r)
  tempR <- read.csv(file.path(results, regions[r], "predictedValues.csv"), header = FALSE)
  
  
  if( !(all(is.na(tempR)))){
    tempR[,1] <- exp(tempR[,1]) - 1
    
    minValues[r] <- min(tempR, na.rm = TRUE)
    maxValues[r] <- max(tempR, na.rm = TRUE)
    
    tempR <- na.omit(tempR)
    
    allValues <- c(allValues, as.vector(tempR[,1]))
  }else{print("All values are NAs")}
  
}

mean(allValues)
median(allValues)
sd(allValues)

quantile(allValues, probs = seq(0.80, 1, by = 0.01))
quantile(allValues, probs = seq(0.99, 1, by = 0.0001))


hist(allValues)
hist(allValues, xlim = c(0, 500), breaks = seq(floor(min(allValues, na.rm = TRUE)),  ceiling(max(allValues, na.rm = TRUE)), by = 1))


minV <- min(minValues, na.rm = TRUE)
maxV <- max(maxValues, na.rm = TRUE)


###########################################

colbrks <-  c(minV, seq(5, 100, length.out = 198), maxV)

r.cols <- magma(199)

png(file.path(figures, "Abundance_corrected.png"),width=17.5,height=8.75,units="cm",res=resdpi)
#pdf(file.path(figures, "Abundance.pdf"),width= wide_inch, height= wide_inch/2, pointsize = point_size)

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
par(mar=c(1,10,1,10))
scale <- c(rep(magma(199)[1], times = 20), rep(magma(199), each = 2), rep(magma(199)[199], times = 20))
b <- barplot(rep(1, 438), col = scale, border =scale, axes = FALSE )
mtext("5", at = b[20], cex = legendcex)
mtext("100", at = b[418], cex = legendcex)
mtext(expression("Abundance (individuals per" ~ m^{2} ~ ")"), side = 1, at = 250, cex = legendcex - 0.2)

dev.off()



