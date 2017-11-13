########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

library(raster)
########################################################
# 1. Load Libraries and Data
########################################################

esa <- raster(file.path("P:\\sWorm_GIS_data\\Land_cover\\esa_land_use_sworm.tif"))
str(esa)
head(esa)

bio10_5 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim\\CHELSA_bio10_5.tif")
bio10_13 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim\\CHELSA_bio10_13.tif")
bio10_14 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim\\CHELSA_bio10_14.tif")

bio10_5_res <- disaggregate(bio10_5, res(bio10_5) / res(esa))
bio10_13_res <- disaggregate(bio10_13, res(bio10_13) / res(esa))
bio10_14_res <- disaggregate(bio10_14, res(bio10_14) / res(esa))

rm(bio10_5); rm(bio10_13); rm(bio10_14)


bio10_5_res <- bio10_5_res/10

# ph1 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\SoilGridsData\\PHIHOX_M_sl1_1km_ll.tif")
# ph2 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\SoilGridsData\\PHIHOX_M_sl2_1km_ll.tif")
ph3 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\SoilGridsData\\PHIHOX_M_sl3_1km_ll.tif")
# ph4 <- raster("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\SoilGridsData\\PHIHOX_M_sl4_1km_ll.tif")

# ph <- stack(ph1, ph2, ph3, ph4)

# ph <- weighted.mean(ph, w = c(1, 5, 15, 30), na.rm = TRUE)  

ph <- disaggregate(ph3, res(ph3) / res(esa))
ph <- ph/10