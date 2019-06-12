## Reviewer suggested to keep the number of sites in each
## band constant
## This script does that


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

sitesbylat <- spp[,c('Study_site', 'Latitude__decimal_degrees')]
sitesbylat <- unique(sitesbylat[c("Study_site", "Latitude__decimal_degrees")])


### How many sites should each band contain? We had 23 bands before...so keep the same?
nrow(sitesbylat) / 23
# 232.7826

## Put in order by latitude
sitesbylat <- sitesbylat[order(sitesbylat$Latitude__decimal_degrees),]


## Ok, we are fine splitting up studies
## What we are not fine today, is splitting up sites 
## that have the same coordinates

let <- 1

sitesbylat$band <- NA
sitesbylat$band[1] <- letters[let]

epsilon <- 20 ## This has worked the best so far!

equalband <- (ceiling(nrow(sitesbylat) / 23)) +20 ## This will mean taht we will be closer to the actual value

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

colNames <- c("band", "Number of Binomials", "Number of Morphospecies", "number of genus", "number of sites", "minLat", "maxLat")
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
  bandDat[i, 5] <- length(unique(bnd$Study_site))
  
  # min lat
  bandDat[i, 'minLat'] <- min(bnd$Latitude__decimal_degrees)
  bandDat[i, 'maxLat'] <- max(bnd$Latitude__decimal_degrees)
  
  
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

## Add in all the gaps




bandDat[23, 1] <- "w"
  
bandDat[23, 2:ncol(bandDat)]<- c(0, 0, 0, 0, -36.999999, -35.399999,  0) # Theres a big gap with no sites

bandDat$latDiff <- bandDat$maxLat - bandDat$minLat

bandDat <- bandDat[order(bandDat$minLat),]


jpeg(file = file.path(figures, "LDG_regional_fixedSites.jpg"), quality = 100, res = 200, width = 2000, height = 1000)

b <- barplot(bandDat$total, width = bandDat$latDiff, space = 0, xaxs = "i", ylab = "Number of Species", xlab = "Latitude")

min(bandDat$minLat)
max(bandDat$maxLat)

## Where would -40 be??
# 0 =  -40.21667
# 107.91 = 68.4525
# so x spans 108.6692 degrees
# each 1  unit on the x = 1.01 degree

#axis(1, at = c(0), labels = c(1))
#axis(1, at = c(107.91), labels = c(2))

## So where is -40 degrees
# -40 degrees is 0.21667 change
1 / (1.01/0.21667)

axis(1, at = 0.214, labels = "-40")

# -30
# 10.21667 degrees change
# 1 / (1.01 / 10.21667)
# axis(1, at = 10.115511, labels = "-30")

## ok, a sequency of every 10 degrees
# from 0.214 to 108, every (10.11551 - 0.214) = 9.90151

# / (1.01/ 8.452500) # 8.45 degree change
#axis(1, at = (99.22910 + 8.368812), labels = "68.4")

labs <- as.character(seq(-40, 60, by = 10))
labs <- c(labs, "68.4")

axis(1, at = c(seq(0.214, 108, by = 9.90151), (99.22910 + 8.368812)), labels = labs)

# axis(1, at = (99.22910 + 9.90151), labels = "70")
dev.off()



