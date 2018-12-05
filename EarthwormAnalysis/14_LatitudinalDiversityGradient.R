########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
library(lme4)
library(car)

source("Functions/FormatData.R")
source(file.path("Functions", "Plots.R"))
# source("Functions/lme4_ModellingFunctions.R")
# source("Functions/ModelSimplification.R")
# source("MEE3_1_sm_Appendix_S1/HighstatLib.R")

#################################################
# 2 Create folders
#################################################

figures <- "Figures"

if(!dir.exists("14_Data")){
  dir.create("14_Data")
}
data_out <- "14_Data"

#################################################
# 3. Loading in data
#################################################
data_in <- "10_Data"
files_spp <- list.files(file.path(data_in))
files_spp <- files_spp[grep("SpecieswithFunctionalGroups_", files_spp)]
file_dates <- sapply(strsplit(files_spp, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files_spp[grep(date, files_spp)]


spp <- read.csv(file.path(data_in, loadin))
rm(files_spp)

data_in <- "4_Data"
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
# 4. Quick and dirty plot
#################################################

plot(sites$Latitude__decimal_degrees, sites$SpeciesRichness)


#################################################
# 5. Model with polynomial
#################################################

sites$polyLatitude <- poly(sites$Latitude__decimal_degrees, 2)

ldg1 <- glm(SpeciesRichness ~  polyLatitude * ExtractionMethod
            , data = sites, family = poisson)
# AIC - 17595.94
ldg2 <- glm(SpeciesRichness ~  polyLatitude + ExtractionMethod
            , data = sites, family = poisson)
# AIC - 17777.55


ldg <- glm(SpeciesRichness ~  poly(Latitude__decimal_degrees, 2) * ExtractionMethod
           , data = sites, family = poisson)
vars <- list()
vars[[1]] <- unique(sites$ExtractionMethod)[2:length(unique(sites$ExtractionMethod))]
vars[[2]] <- seq(min(sites$Latitude__decimal_degrees), max(sites$Latitude__decimal_degrees), by = 0.1)
newdata <- expand.grid(vars)
names(newdata) <- c("ExtractionMethod", "Latitude__decimal_degrees")
newdata$predicted <- predict(ldg, newdata,type="response")

handsorting <- newdata[newdata$ExtractionMethod == "Hand sorting",]
handsorting <- handsorting[handsorting$Latitude__decimal_degrees > min(sites$Latitude__decimal_degrees[sites$ExtractionMethod == "Hand sorting"], na.rm = TRUE),]
handsorting <- handsorting[handsorting$Latitude__decimal_degrees < max(sites$Latitude__decimal_degrees[sites$ExtractionMethod == "Hand sorting"], na.rm = TRUE),]

mustard <- newdata[newdata$ExtractionMethod == "Liquid extraction (e.g. Mustard)",]
mustard <- mustard[mustard$Latitude__decimal_degrees > min(sites$Latitude__decimal_degrees[sites$ExtractionMethod == "Liquid extraction (e.g. Mustard)"], na.rm = TRUE),]
mustard <- mustard[mustard$Latitude__decimal_degrees < max(sites$Latitude__decimal_degrees[sites$ExtractionMethod == "Liquid extraction (e.g. Mustard)"], na.rm = TRUE),]

handsortingandmustard <- newdata[newdata$ExtractionMethod == "Hand sorting + Liquid extraction (e.g. Mustard)",]
handsortingandmustard <- handsortingandmustard[handsortingandmustard$Latitude__decimal_degrees > min(sites$Latitude__decimal_degrees[sites$ExtractionMethod == "Hand sorting + Liquid extraction (e.g. Mustard)"], na.rm = TRUE),]
handsortingandmustard <- handsortingandmustard[handsortingandmustard$Latitude__decimal_degrees < max(sites$Latitude__decimal_degrees[sites$ExtractionMethod == "Hand sorting + Liquid extraction (e.g. Mustard)"], na.rm = TRUE),]

jpeg(file = file.path(figures, "LDG_local.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
plot(handsorting$Latitude__decimal_degrees, handsorting$predicted, pch = 16, xlab = "Latitude", ylab = "Species Richness", type = "l", ylim = c(0,  3.5))
lines(mustard$Latitude__decimal_degrees, mustard$predicted, col = "red")
lines(handsortingandmustard$Latitude__decimal_degrees, handsortingandmustard$predicted, col = "blue")
legend("topleft", legend= c("Handsorting", "Mustard", "Both"), col = c("black", "red", "blue"), lwd = 1, bty = "n")
dev.off()
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
    bandDat[i, 2] <- length(unique(bnd$Revised))
    
    
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


RichnessDat <- bandDat
###############################################
## Trying to correct for area
################################################
library(rgdal)

area <- readOGR(dsn = "I:/sWorm/Database/Area", layer = "area_5deg_equal_terrestrial_zero")
dat <- area@data

sum(dat$Shape_Area)
dat$prop <- dat$Shape_Area / sum(dat$Shape_Area)

dat <- dat[-(1:9) ,]
dat <- dat[dat$Id %in% 1:32,]

bandDat <- cbind(bandDat, dat)

#### A linear model
# What is the effect of latitude after accounting for 1) area and 2) number of sampling sites

lm1 <- glm(total ~ Shape_Area, family = poisson, data = bandDat)
bandDat$area_Resid <- resid(lm1)
plot(y = bandDat$area_Resid, x = 1:nrow(bandDat))
axis(1,  labels = bandDat$band, at = 1:nrow(bandDat))

lm2 <- glm(total ~ bandDat$`number of sites`, family = poisson, data = bandDat)
bandDat$nSites_Resid <- resid(lm2)

plot(y = bandDat$nSites_Resid, x = 1:nrow(bandDat))
axis(1,  labels = bandDat$band, at = 1:nrow(bandDat))


#########################
## Calculating rarefied richness
############################
sitesbyband <- spp[,c('Study_site', 'band')]
sitesbyband <- unique(sitesbyband[c("Study_site", "band")])

# n_min <- min(table(sitesbyband$band)[(which(table(sitesbyband$band) != 0))]) ## the latitude with the fewest number of samples/sites
n_min <- 40 # number of sites in a latitude, that isn't ridiculously low
n_rarefied <- 1000


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
  if(length(unique(bnd$Study_site)) < n_min){
    
    bandDat[i, 3] <- 0
    }else{
    
    ## These are the unique sites in the band. 
    ## Take names of random sample to split the spp dataframe by
    for(r in 1:n_rarefied){
      sbyb <- sitesbyband[sitesbyband$band == levels(spp$band)[i],]
      thesesites <- sample(sbyb$Study_site, n_min)
      
      bnd_temp <- bnd[bnd$Study_site %in% thesesites,]
      
      
      # Number of species binomials
      rarefiedDat[r, 1] <- length(unique(bnd_temp$Revised))
      
      
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

par(mar = c(4, 4, 1, 5))
barplot(bandDat$rarefiedRichness, space = 0, xaxs = "i", ylab = "Rarefied Richness", xlab = "Latitude")
axis(1, at = 0:23, labels = seq(-45, 70, by = 5))
dev.off()     

  