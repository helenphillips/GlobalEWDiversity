########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

#################################################
# 1. Loading libraries
#################################################
library(sp)
library(raster)

# library(maptools)
# library(maps)
# library(plyr)
# library(dplyr)

#################################################
# 2. Loading in variables
#################################################

data_in <-"0_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinsites <- loadin[grep("sites_", loadin)]

if(!dir.exists("1_Data")){
  dir.create("1_Data")
}
data_out <- "1_Data"

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

#################################################
# 3. Load in data
#################################################

sites <- read.csv(file.path(data_in, loadinsites))
# sites <- SiteLevels(sites) 

#################################################
# 4. Load in TIFs
#################################################
## This would need to loop through all the diffferent tifs

tifs <- c("bio10_1","bio10_2","bio10_3","bio10_4", "bio10_5", 
          "bio10_6", "bio10_7","bio10_8","bio10_9","bio10_10",
          "bio10_11","bio10_12", "bio10_13", "bio10_14","bio10_15",
          "bio10_16","bio10_17","bio10_18","bio10_19")
divide <- tifs[c(1, 5:11)]
for(t in tifs){
  tif <- raster(file.path("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim", 
                  paste("CHELSA_", t, ".tif", sep="")))
  sites$V1 <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees , sites$Latitude__decimal_degrees))
    if(t %in% divide){
      sites$V1 <- sites$V1 / 10
    }
  names(sites)[names(sites) == "V1"] <- t
}

#### 
## Check that any without values are because they have no coordinate
#################################################
# 5. Save data
#################################################

write.csv(sites, file = file.path(data_out, paste("sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
