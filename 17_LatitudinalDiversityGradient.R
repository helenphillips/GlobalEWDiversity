########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
library(lme4)
library(car)

source(file.path("Functions", "FormatData.R"))
source(file.path("Functions", "Plots.R"))
# source("Functions/lme4_ModellingFunctions.R")
# source("Functions/ModelSimplification.R")
# source("MEE3_1_sm_Appendix_S1/HighstatLib.R")

#################################################
# 2 Create folders
#################################################

figures <- "Figures"

if(!dir.exists("17_Data")){
  dir.create("17_Data")
}
data_out <- "17_Data"

#################################################
# 3. Loading in data
#################################################
data_in <- "16_Data"
files_spp <- list.files(file.path(data_in))
files_spp <- files_spp[grep("SpecieswithFunctionalGroups_", files_spp)]
file_dates <- sapply(strsplit(files_spp, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files_spp[grep(date, files_spp)]


spp <- read.csv(file.path(data_in, loadin))
rm(files_spp)

data_in <- "8_Data"
files_sites <- list.files(file.path(data_in))
files_sites <- files_sites[grep("sitesRichness_", files_sites)]
file_dates <- sapply(strsplit(files_sites, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files_sites[grep(date, files_sites)]

sites <- read.csv(file.path(data_in, loadin))

rm(files_sites)
rm(date)
rm(loadin)

#################################################
# 6. LDG by bands
#################################################
# summary(spp$Latitude__decimal_degrees)

# Want the sites that we used in the modelling
spp <- spp[spp$Study_site %in% unique(sites$Study_site),]


brks <- seq(-45, 70, by = 5)
spp$band <- cut(spp$Latitude__decimal_degrees, breaks = brks, labels = NULL)

bandDat <- data.frame(matrix(rep(NA, times = 5 * length(levels(spp$band))), 
                             ncol = 5, nrow  = length(levels(spp$band))))

names(bandDat) <- c("band", "Number of Binomials", "Number of Morphospecies", "number of genus", "number of sites")

for(i in 1:length(levels(spp$band))){
  bnd <- spp[spp$band == levels(spp$band)[i],] 
  print(levels(spp$band)[i])
  
  # Which band
  bandDat[i, 1] <- levels(spp$band)[i]
  # Sampling effort
  bandDat[i, 5] <- length(unique(bnd$Study_site))
  
  if(nrow(bnd) == 0){
    bandDat[i, c(2:4)] <- 0
  }else{
    
    # Number of species binomials
    
    
    binomials <- unique(bnd$Revised)
    binomials <- binomials[!(is.na(binomials))]
    
    
    bandDat[i, 2] <- length(binomials)
    
    
    # Number of morphospecies
    mrphs <- bnd[!(is.na(bnd$MorphospeciesID)),]
    if(nrow(mrphs) > 0){
      bandDat[i, 3] <- length(unique(paste(mrphs$Genus, mrphs$MorphospeciesID)))
    }else{ bandDat[i, 3] <- 0}
    
    # Number of genus not considered anywhere else
    gns <- bnd[is.na(bnd$SpeciesBinomial),] 
    gns <- gns[is.na(gns$MorphospeciesID),] 
    uni <- as.character(unique(gns$Genus))
    if(length(uni) > 0){
      for(g in 1:length(uni)){
        matches <- grep(uni[g], bnd$Revised)
        if(length(matches) > 0){
          uni[g] <- 0
        } else {uni[g] <- 1}
      }
      uni <- as.numeric(uni)
      bandDat[i, 4] <- sum(uni)
    } else {bandDat[i, 4] <- 0}
  }
}


bandDat$total <- rowSums(bandDat[,2:4])

jpeg(file = file.path(figures, "LDG_regional.jpg"), quality = 100, res = 200, width = 2000, height = 1000)

par(mar = c(4, 4, 1, 5))
b <- barplot(bandDat$total, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")
barplot(bandDat$total, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")
par(new=TRUE)
plot(b[,1],bandDat[,5],xaxs = "i", xlim=c(0,23),type="l",col="red",axes=FALSE,ylim=c(0,1200),ann=FALSE)
axis(4,at=seq(0,1200,100), las = 2)
axis(1, at = 0:23, labels = seq(-45, 70, by = 5))
mtext("Number of sites", side = 4, line = 3)
dev.off()     

### Without morphospecies

bandDat$totalnoMorphs <- rowSums(bandDat[,c(2,4)])
jpeg(file = file.path(figures, "LDG_regional_nomorphs.jpg"), quality = 100, res = 200, width = 2000, height = 1000)

par(mar = c(4, 4, 1, 5))
b <- barplot(bandDat$totalnoMorphs, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")
barplot(bandDat$totalnoMorphs, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")
par(new=TRUE)
plot(b[,1],bandDat[,5],xaxs = "i", xlim=c(0,23),type="l",col="red",axes=FALSE,ylim=c(0,1200),ann=FALSE)
axis(4,at=seq(0,1200,100), las = 2)
axis(1, at = 0:23, labels = seq(-45, 70, by = 5))
mtext("Number of sites", side = 4, line = 3)
dev.off()     



RichnessDat <- bandDat


#########################
## Calculating rarefied richness
############################
sitesbyband <- spp[,c('Study_site', 'band')]
sitesbyband <- unique(sitesbyband[c("Study_site", "band")])

# n_min <- min(table(sitesbyband$band)[(which(table(sitesbyband$band) != 0))]) ## the latitude with the fewest number of samples/sites
n_min <- 22 # number of sites in a latitude, that isn't ridiculously low
n_rarefied <- 1000


bandDat <- data.frame(matrix(rep(NA, times = 4 * length(levels(spp$band))), 
                             ncol = 4, nrow  = length(levels(spp$band))))
names(bandDat) <- c("band", "number of sites", "rarefiedRichness", 'rarefiedNoMorphs')


rarefiedDat <- data.frame(matrix(rep(NA, times = 5 * n_rarefied), 
                                 ncol = 5, nrow  = n_rarefied))
names(rarefiedDat) <- c("Number of Binomials", "Number of Morphospecies", "number of genus", "total", "totalnoMorphs")

for(i in 1:length(levels(spp$band))){
  bnd <- spp[spp$band == levels(spp$band)[i],] 
  print(levels(spp$band)[i])
  
  # Which band
  bandDat[i, 1] <- levels(spp$band)[i]
  # Sampling effort
  bandDat[i, 2] <- length(unique(bnd$Study_site))
  
  if(nrow(bnd) == 0){
    bandDat[i, 3] <- 0
    bandDat[i, 4] <- 0    }
  if(length(unique(bnd$Study_site)) < n_min){
    
    bandDat[i, 3] <- 0
    bandDat[i, 4] <- 0
    }else{
    
    ## These are the unique sites in the band. 
    ## Take names of random sample to split the spp dataframe by
    for(r in 1:n_rarefied){
      sbyb <- sitesbyband[sitesbyband$band == levels(spp$band)[i],]
      thesesites <- sample(sbyb$Study_site, n_min)
      
      bnd_temp <- bnd[bnd$Study_site %in% thesesites,]
      
      
      # Number of species binomials
      binomials <- unique(bnd_temp$Revised)
      binomials <- binomials[!(is.na(binomials))]
      
      
      rarefiedDat[r, 1] <- length(binomials)
      
      
      # Number of morphospecies
      mrphs <- bnd_temp[!(is.na(bnd_temp$MorphospeciesID)),]
      if(nrow(mrphs) > 0){
        rarefiedDat[r, 2] <- length(unique(paste(mrphs$Genus, mrphs$MorphospeciesID)))
      }else{ rarefiedDat[r, 2] <- 0}
      
      # Number of genus not considered anywhere else
      gns <- bnd_temp[is.na(bnd_temp$SpeciesBinomial),] 
      gns <- gns[is.na(gns$MorphospeciesID),] 
      uni <- as.character(unique(gns$Genus))
      if(length(uni) > 0){
        for(g in 1:length(uni)){
          matches <- grep(uni[g], bnd_temp$Revised)
          if(length(matches) > 0){
            uni[g] <- 0
          } else {uni[g] <- 1}
        }
        uni <- as.numeric(uni)
        rarefiedDat[r, 3] <- sum(uni)
      } else {rarefiedDat[r, 3] <- 0}
      
      # print(mean(rarefiedDat$`Number of Morphospecies`))
      
      rarefiedDat$total <- rowSums(rarefiedDat[,1:3])
      rarefiedDat$totalnoMorphs <- rowSums(rarefiedDat[,c(1,3)])
    }
    ## Done the 1000 samples
    bandDat[i, 3] <- mean(rarefiedDat$total)
    bandDat[i, 4] <- mean(rarefiedDat$totalnoMorphs)
  }
}

jpeg(file = file.path(figures, "LDG_regional_paired.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
par(mfrow = c(1, 2))
par(mar = c(4, 4, 1, 5))
b <- barplot(RichnessDat$total, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")
# barplot(RichnessDat$total, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")
par(new=TRUE)
plot(b[,1],RichnessDat[,5],xaxs = "i", xlim=c(0,23),type="l",col="red",axes=FALSE,ylim=c(0,1200),ann=FALSE)
axis(4,at=seq(0,1200,100), las = 2)
axis(1, at = 0:23, labels = seq(-45, 70, by = 5))
mtext("Number of sites", side = 4, line = 3)
mtext("(a)", side = 3, line = 0, at = 0, adj = 0.1)

par(mar = c(4, 4, 1, 5))
barplot(bandDat$rarefiedRichness, space = 0, xaxs = "i", ylab = "Rarefied Richness", xlab = "Latitude")
axis(1, at = 0:23, labels = seq(-45, 70, by = 5))
mtext("(b)", side = 3, line = 0, at = 0, adj = 0.1)

dev.off()     

jpeg(file = file.path(figures, "LDG_regional_paired_noMorphs.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
par(mfrow = c(1, 2))
par(mar = c(4, 4, 1, 5))
b <- barplot(RichnessDat$totalnoMorphs, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")
# barplot(RichnessDat$totalnoMorphs, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")
par(new=TRUE)
plot(b[,1],RichnessDat[,5],xaxs = "i", xlim=c(0,23),type="l",col="red",axes=FALSE,ylim=c(0,1200),ann=FALSE)
axis(4,at=seq(0,1200,100), las = 2)
axis(1, at = 0:23, labels = seq(-45, 70, by = 5))
mtext("Number of sites", side = 4, line = 3)
mtext("(a)", side = 3, line = 0, at = 0, adj = 0.1)

par(mar = c(4, 4, 1, 5))
barplot(bandDat$rarefiedNoMorphs, space = 0, xaxs = "i", ylab = "Rarefied Richness (excluding morphospecies)", xlab = "Latitude")
axis(1, at = 0:23, labels = seq(-45, 70, by = 5))
mtext("(b)", side = 3, line = 0, at = 0, adj = 0.1)

dev.off()     


##################################################################################
# Doing the whole things for different n_mins
##################################################################################

n_min <- c(3, 5, 9, 22, 26, 30, 40, 47, 55) # number of sites in a latitude, that isn't ridiculously low
n_rarefied <- 1000

jpeg(file = file.path(figures, "LDG_regional_allnmins.jpg"), quality = 100, res = 200, width = 4000, height = 4000)
par(mfrow = c(3, 3))

for(m in n_min){

bandDat <- data.frame(matrix(rep(NA, times = 3 * length(levels(spp$band))), 
                             ncol = 3, nrow  = length(levels(spp$band))))
names(bandDat) <- c("band", "number of sites", "rarefiedRichness")


rarefiedDat <- data.frame(matrix(rep(NA, times = 4 * n_rarefied), 
                                 ncol = 4, nrow  = n_rarefied))
names(rarefiedDat) <- c("Number of Binomials", "Number of Morphospecies", "number of genus", "total")

for(i in 1:length(levels(spp$band))){
  bnd <- spp[spp$band == levels(spp$band)[i],] 
  print(levels(spp$band)[i])
  
  # Which band
  bandDat[i, 1] <- levels(spp$band)[i]
  # Sampling effort
  bandDat[i, 2] <- length(unique(bnd$Study_site))
  
  if(nrow(bnd) == 0){
    bandDat[i, 3] <- 0}
  if(length(unique(bnd$Study_site)) < m){
    
    bandDat[i, 3] <- 0
  }else{
    
    ## These are the unique sites in the band. 
    ## Take names of random sample to split the spp dataframe by
    for(r in 1:n_rarefied){
      sbyb <- sitesbyband[sitesbyband$band == levels(spp$band)[i],]
      thesesites <- sample(sbyb$Study_site, m)
      
      bnd_temp <- bnd[bnd$Study_site %in% thesesites,]
      
      
      # Number of species binomials
     
      binomials <- unique(bnd_temp$Revised)
      binomials <- binomials[!(is.na(binomials))]
      
      
      rarefiedDat[r, 1] <- length(binomials)
      
      # Number of morphospecies
      mrphs <- bnd_temp[!(is.na(bnd_temp$MorphospeciesID)),]
      if(nrow(mrphs) > 0){
        rarefiedDat[r, 2] <- length(unique(paste(mrphs$Genus, mrphs$MorphospeciesID)))
      }else{ rarefiedDat[r, 2] <- 0}
      
      # Number of genus not considered anywhere else
      gns <- bnd_temp[is.na(bnd_temp$SpeciesBinomial),] 
      gns <- gns[is.na(gns$MorphospeciesID),] 
      uni <- as.character(unique(gns$Genus))
      if(length(uni) > 0){
        for(g in 1:length(uni)){
          matches <- grep(uni[g], bnd_temp$Revised)
          if(length(matches) > 0){
            uni[g] <- 0
          } else {uni[g] <- 1}
        }
        uni <- as.numeric(uni)
        rarefiedDat[r, 3] <- sum(uni)
      } else {rarefiedDat[r, 3] <- 0}
      
      rarefiedDat$total <- rowSums(rarefiedDat[,1:3])
    }
    ## Done the 1000 samples
    bandDat[i, 3] <- mean(rarefiedDat$total)
  }
}



par(mar = c(4, 4.5, 3.3, 5))
barplot(bandDat$rarefiedRichness, space = 0, xaxs = "i", 
        ylab = paste("Rarefied Richness (n_min = ", m,")", sep = ""), xlab = "Latitude")
axis(1, at = 0:23, labels = seq(-45, 70, by = 5))
mtext(paste("n = ", m, sep = ""), side = 3, line = 0, at = 0, adj = 0.1, cex = 3)

}
dev.off()     
