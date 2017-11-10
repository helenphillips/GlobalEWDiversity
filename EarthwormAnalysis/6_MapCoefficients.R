########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}


########################################################
# 1. Load Libraries and Data
########################################################

esa <- raster(file.path("P:\\sWorm_GIS_data\\Land_cover\\esa_land_use_sworm.tif"))
str(esa)
head(esa)

bio10_5 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim\\CHELSA_bio10_5.tif")
bio10_13 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim\\CHELSA_bio10_13.tif")
bio10_14 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim\\CHELSA_bio10_14.tif")

bio10_5_res <- disaggregate(bio10_5, floor(res(bio10_5) / res(esa)))
