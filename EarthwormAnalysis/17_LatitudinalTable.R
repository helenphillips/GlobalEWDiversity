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

Colnames <- c("band", "number of sites", "Number of Binomials", "meanLatRange",  "Number of Morphospecies", "percent native", "percent non-native", "percent unknown")


bandDat <- data.frame(matrix(rep(NA, times = length(Colnames) * length(levels(spp$band))), 
                             ncol = length(Colnames), nrow  = length(levels(spp$band))))

names(bandDat) <- Colnames

for(i in 1:length(levels(spp$band))){
  bnd <- spp[spp$band == (levels(spp$band))[i],] 
  print(levels(spp$band)[i])
  
  # Which band
  bandDat[i, 'band'] <- (levels(spp$band))[i]
  # Sampling effort
  bandDat[i, 'number of sites'] <- length(unique(bnd$Study_site))
  
  if(nrow(bnd) == 0){
    # bandDat[i, c(2:4)] <- 0
    next
  }else{
    
    # Number of species binomials
    
    binomials <- unique(bnd$Revised)
    binomials <- binomials[!(is.na(binomials))]
    
    
    bandDat[i, 'Number of Binomials'] <- length(binomials)
    
    spp$Revised <- as.character(spp$Revised)
    binomials <- as.character(binomials)
    
    ranges <- c()
    if(length(binomials) > 0){
      
      for(b in 1:length(binomials)){
     
        print(b)
      
        those_species <- droplevels(spp[which(spp$Revised == binomials[b]),])
        min_lat <- min(those_species$Latitude__decimal_degrees, na.rm = TRUE)
        print(min_lat)
        max_lat <- max(those_species$Latitude__decimal_degrees, na.rm = TRUE)
        print(max_lat)
      
        range <-  max_lat - min_lat
        ranges <- c(ranges, range)
      }
      
      bandDat[i, 'meanLatRange'] <- mean(ranges)
      
    }
    
    
    
    
    
    # Number of morphospecies
    mrphs <- bnd[!(is.na(bnd$MorphospeciesID)),]
    if(nrow(mrphs) > 0){
      bandDat[i, 'Number of Morphospecies'] <- length(unique(paste(mrphs$Genus, mrphs$MorphospeciesID)))
    }else{ bandDat[i, 'Number of Morphospecies'] <- 0}
    
    # Number of genus not considered anywhere else
    # gns <- bnd[is.na(bnd$SpeciesBinomial),] 
    # gns <- gns[is.na(gns$MorphospeciesID),] 
    # uni <- as.character(unique(gns$Genus))
    # if(length(uni) > 0){
    #   for(g in 1:length(uni)){
    #     matches <- grep(uni[g], bnd$Revised)
    #     if(length(matches) > 0){
    #       uni[g] <- 0
    #     } else {uni[g] <- 1}
    #   }
    #   uni <- as.numeric(uni)
    #   bandDat[i, 4] <- sum(uni)
    # } else {bandDat[i, 4] <- 0}
    
    print(table(bnd$Native.Nonnative))
    bandDat[i, 'percent native'] <- as.vector((table(bnd$Native.Nonnative)[1] / nrow(bnd)) * 100)
    bandDat[i, 'percent non-native'] <- as.vector((table(bnd$Native.Nonnative)[2] / nrow(bnd)) * 100)
    n_unknown <- nrow(bnd) - as.vector((table(bnd$Native.Nonnative)[1] + table(bnd$Native.Nonnative)[2] ))
    
    bandDat[i, 'percent unknown'] <- (n_unknown / nrow(bnd) * 100)

  }
}

bandDat <- bandDat[nrow(bandDat):1, ]

bandDat$meanLatRange <- round(bandDat$meanLatRange, 2)
bandDat[,'percent native'] <- round(bandDat[,'percent native'], 2)
bandDat$'percent non-native' <- round(bandDat$'percent non-native', 2)
bandDat$'percent unknown' <- round(bandDat$'percent unknown', 2)

write.csv(bandDat, file = file.path(data_out, "LatitudinalTable.csv"), row.names = FALSE)
