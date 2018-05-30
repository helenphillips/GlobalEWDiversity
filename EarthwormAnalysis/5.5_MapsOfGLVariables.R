########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

########################################################
# 1. Libraries
########################################################


library(rgdal)
library(raster)
library(rasterVis)

########################################################
# 2. Load In
########################################################
# Habitat Cover
esa <- raster('I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\ESA HC\\LC_1000_NEW\\LC_1000_NEW\\ESACCI-LC-L4-LCCS-Map-1000m-P5Y-2010-v1.6.1_NEW.tif')
plot(esa)
# This would need a legend


## Soil
soil <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\SoilGrids\\1km"
pH <- raster(file.path(soil, "PHIHOX_weighted.tif"))
plot(pH)

par(mar=c(5, 4, 1, 1))
levelplot(pH, par.settings = RdBuTheme, margin=FALSE)



OrgC <- raster(file.path(soil, "ORCDRC_weighted.tif"))
plot(OrgC)


clay <- raster(file.path(soil, "CLYPPT_weighted.tif"))
plot(clay)

silt <- raster(file.path(soil, "SLTPPT_weighted.tif")) 
plot(silt)

CEC <- raster(file.path(soil, "CECSOL_weighted.tif"))
plot(CEC)

## Climate
chelsa <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\CHELSAData\\BioClim"

AnnPrecip <- raster(file.path(chelsa, "CHELSA_bio10_1.tif"))
plot(AnnPrecip)


PrecipSeason <- raster(file.path(chelsa, "CHELSA_bio10_15.tif"))
plot(PrecipSeason)


snow <- raster("I:\\sWorm\\Database\\snow\\snow_2015_sum.tif")
plot(snow)

PET <- raster("I:\\sWorm\\Database\\PET_he_annual\\pet_he_yr_TIF.tif")
plot(PET)

PETSD <- raster("I:\\sWorm\\Database\\PET_he_monthly\\pet_he_SD.tif")
plot(PETSD)

AI <- raster("I:\\sWorm\\Database\\AI_annual\\ai_yr_TIF.tif")
plot(AI)