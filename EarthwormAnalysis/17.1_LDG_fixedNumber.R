## Reviewer suggested to keep the number of sites in each
## band constant
## This script does that


########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "TSGIS02"){
  setwd("C:/sWorm/EarthwormAnalysis")
}

#######################################
# variables
######################################

wide_cm <- 12
wide_inch <- 4.75
wide_inch_small <- 2.25
point_size <- 7
plotlabcex <- 0.5
resdpi <- 300


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

sitesbylat <- spp[,c('Study_site', 'Latitude__decimal_degrees')]
sitesbylat <- unique(sitesbylat[c("Study_site", "Latitude__decimal_degrees")])


### How many sites should each band contain? We had 23 bands before...so keep the same?
nrow(sitesbylat) / 23
# 232.8696

## Put in order by latitude
sitesbylat <- sitesbylat[order(sitesbylat$Latitude__decimal_degrees),]


## Ok, we are fine splitting up studies
## What we are not fine with, is splitting up sites 
## that have the same coordinates

let <- 1

sitesbylat$band <- NA
sitesbylat$band[1] <- letters[let]

epsilon <- 20 ## This has worked the best so far!

equalband <- (ceiling(nrow(sitesbylat) / 23)) + 25 ## This will mean that we will be closer to the actual value

for(r in 2:nrow(sitesbylat)){
  ## If there's not many in the band, then we can assign the same letter
  if(nrow(sitesbylat[which(sitesbylat$band == letters[let]),]) < ( equalband- epsilon )){
    sitesbylat$band[r] <- letters[let]
    next
  }else{
    if(sitesbylat$Latitude__decimal_degrees[r] == sitesbylat$Latitude__decimal_degrees[r-1]){
      sitesbylat$band[r] <- letters[let]
      next
    }
    else{let <- let + 1
    sitesbylat$band[r] <- letters[let]
    next} # if they are not the same, then advance the letter
    
  }
}
  

sitesbylat$Latitude__decimal_degrees <- NULL
spp <- merge(spp, sitesbylat, by = "Study_site", all.x = TRUE)

colNames <- c("band", "Number of Binomials", "Number of Morphospecies", "number of genus", 
              "number of genus no morphs", "number of sites", "minLat", "maxLat")
bandDat <- data.frame(matrix(rep(NA, times = length(colNames) * length(levels(spp$band))), 
                             ncol = length(colNames), nrow  = length(levels(spp$band))))

names(bandDat) <- colNames




spp$band <- as.factor(spp$band)

for(i in 1:length(levels(spp$band))){
  bnd <- spp[spp$band == levels(spp$band)[i],] 
  print(levels(spp$band)[i])
  
  # Which band
  bandDat[i, 1] <- levels(spp$band)[i]
  # Sampling effort
  bandDat[i, "number of sites"] <- length(unique(bnd$Study_site))
  
  # min lat
  bandDat[i, 'minLat'] <- min(bnd$Latitude__decimal_degrees)
  bandDat[i, 'maxLat'] <- max(bnd$Latitude__decimal_degrees)
  
  
  if(nrow(bnd) == 0){
    bandDat[i, c(2:5)] <- 0
  }else{
    
    # Number of species binomials
    binomials <- unique(bnd$Revised)
    binomials <- binomials[!(is.na(binomials))]
    
    
    bandDat[i, 2] <- length(binomials)
    
    # Number of morphospecies
    mrphs <- bnd[!(is.na(bnd$MorphospeciesID)),]
    if(nrow(mrphs) > 0){
      bandDat[i, "Number of Morphospecies"] <- length(unique(paste(mrphs$Genus, mrphs$MorphospeciesID)))
    }else{ bandDat[i, "Number of Morphospecies"] <- 0}
    
    # Number of genus not considered anywhere else
    gns <- bnd[is.na(bnd$Revised),] 
    gns <- gns[is.na(gns$MorphospeciesID),] 
    uni <- unique(gns$Genus)
    uni <- as.character(uni[!(is.na(uni))])
    
    if(length(uni) > 0){
      for(g in 1:length(uni)){
        matches <- grep(uni[g], bnd$Revised)
        if(length(matches) > 0){
          uni[g] <- 0
        } else {uni[g] <- 1}
      }
      uni <- as.numeric(uni)
      bandDat[i, "number of genus"] <- sum(uni)
    } else {bandDat[i, "number of genus"] <- 0}
    
    # Number of genus not considered else where, if we forget about morphospecies
    
    # Number of genus not considered anywhere else
    gns <- bnd[is.na(bnd$Revised),] 
    uni <- unique(gns$Genus)
    uni <- as.character(uni[!(is.na(uni))])
    
    if(length(uni) > 0){
      for(g in 1:length(uni)){
        matches <- grep(uni[g], bnd$Revised)
        if(length(matches) > 0){
          uni[g] <- 0
        } else {uni[g] <- 1}
      }
      uni <- as.numeric(uni)
      bandDat[i, "number of genus no morphs"] <- sum(uni)
    } else {bandDat[i, "number of genus no morphs"] <- 0}
    
    
  }
}


bandDat$total <- rowSums(bandDat[,c("Number of Binomials", "Number of Morphospecies", "number of genus")])
bandDat$total_no_morphs <- rowSums(bandDat[,c("Number of Binomials", "number of genus no morphs")])

## Add in all the gaps




bandDat[23, 1] <- "w"
  
bandDat[23, 2:ncol(bandDat)]<- c(0, 0, 0, 0, 0, -36.999999, -35.399999, 0, 0) # Theres a big gap with no sites

bandDat$latDiff <- bandDat$maxLat - bandDat$minLat

bandDat <- bandDat[order(bandDat$minLat),]


# jpeg(file = file.path(figures, "LDG_regional_fixedSites.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
pdf(file.path(figures, "LDG_regional_fixedSites.pdf"),width= wide_inch_small, height= wide_inch_small/2, pointsize = point_size)

b <- barplot(bandDat$total, width = bandDat$latDiff, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")

min(bandDat$minLat)
max(bandDat$maxLat)

## Where would -40 be??
# 0 =  -40.21667
# 107 = 68.4525
# so x spans 108.64 degrees
# each 1  unit on the x = 1.015 degree
oneunit <- (max(bandDat$maxLat) - min(bandDat$minLat)) / 107
#axis(1, at = c(0), labels = c(1))
# axis(1, at = c(107), labels = c(2))

## So where is -40 degrees
# -40 degrees is 0.21667 change
minus40 <- 1 / (oneunit/0.21667)

axis(1, at = minus40, labels = "-40")

# -30
# 10.21667 degrees change
# 1 / (oneunit / 10.21667)
# axis(1, at = 10.16519, labels = "-30")

## ok, a sequency of every 10 degrees
# from 0.214 to 107, every (10.16519 - minus40) = 9.949612

# 1 / (oneunit/ 8.43075) # 8.388272 degree change
#axis(1, at = (99.949612 + 8.388272), labels = "68.43")



labs <- as.character(seq(-40, 60, by = 10))
labs <- c(labs, "68.4")

axis(1, at = c(seq(minus40, 107, by = ( 10.05974 - minus40) ), 107), labels = labs)


# axis(1, at = (99.22910 + 9.90151), labels = "70")
dev.off()

#######################
# And with no morphospecies
#########################

jpeg(file = file.path(figures, "LDG_regional_fixedSites_nomorphs.jpg"), quality = 100, res = 200, width = 2000, height = 1000)

b <- barplot(bandDat$total_no_morphs, width = bandDat$latDiff, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")
min(bandDat$minLat)
max(bandDat$maxLat)
oneunit <- (max(bandDat$maxLat) - min(bandDat$minLat)) / 107
minus40 <- 1 / (oneunit/0.21667)

axis(1, at = minus40, labels = "-40")

labs <- as.character(seq(-40, 60, by = 10))
labs <- c(labs, "68.4")

axis(1, at = c(seq(minus40, 107, by = ( 10.05974 - minus40) ), 107), labels = labs)

dev.off()

