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
library(rgdal)

#################################################
# 2. Loading in variables
#################################################

data_in <-"1_Data"

files <- list.files(file.path(data_in))
files <- files[grep("sitesWithChelsaAndSoil_", files)]

file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinsites <- loadin[grep("sitesWithChelsaAndSoil_", loadin)]

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

############## SNOW

snow <- "I:\\sWorm\\Database\\snow"
tif <- raster(file.path(snow, "snow_2015_sum.tif"))
sites$SnowMonths <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees))


############## ARIDITY
arid <- "I:\\sWorm\\Database\\AI_annual"
tif <- raster(file.path(arid, "ai_yr_TIF.tif"))
sites$Aridity <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees))

############## PET

pet <- "I:\\sWorm\\Database\\PET_he_annual"
tif <- raster(file.path(pet, "pet_he_yr_TIF.tif"))
sites$PETyr <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees))

petSD

#################################################
# 5. Save data
#################################################

write.csv(sites, file = file.path(data_out, paste("sitesWithChelsaAndSoilAndOthers_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

