########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\USers\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
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
files <- files[grep("sitesWithChelsa_", files)]

file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinsites <- loadin[grep("sitesWithChelsa_", loadin)]

if(!dir.exists("2_Data")){
  dir.create("2_Data")
}
data_out <- "2_Data"

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


soil <- "I:\\sDiv\\Phillips\\sWorm\\SpatialAnalysis\\SoilGrids\\1km"

tifs <- c("PHIHOX", "CLYPPT", "SLTPPT", "SNDPPT", "CECSOL", "ORCDRC", "TAXNWRB_1")
layers <- c("sl1","sl2", "sl3", "sl4")

filenames <- list.files(soil)

for(t in tifs){
  dl <- filenames[grep(t, filenames)]
  if(t != "TAXNWRB_1"){
    dl <- dl[grep(paste(layers,collapse="|"), dl)] ## Just the layers for the first 30cm
  }
  
  temp <- as.data.frame(rep(NA, length=nrow(sites)))
  
  for(f in 1:length(dl)){
    tif <- raster(file.path(soil, dl[f]))
    
    temp[,f] <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees , sites$Latitude__decimal_degrees))
    
  }
  ##  The soil taxonomy is only one layer
  if(t != "TAXNWRB_1"){
    weight <- c(0.001, 0.05, 0.1, 0.15)
    sites$V1 <- apply(temp[1:4], 1, weighted.mean, weight)
    
    ## PH is multiplied by 10, organic carbon is currently grams per kg
    if(t == "PHIHOX" | t == "ORCDRC"){
      sites$V1 <- sites$V1 / 10
    }
  }else{sites$V1 <- temp[,1]}
  names(sites)[names(sites) == "V1"] <- t
}


# soil_tax <- read.delim(file = "C:\\Users\\Dropbox\\sWorm\\SoilGrids\\1km\\TAXNWRB.txt")
## This has not been sorted

#################################################
# 5. Save data
#################################################

write.csv(sites, file = file.path(data_out, paste("sitesWithChelsaAndSoil_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

