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



if(!dir.exists("Prelim_Figures")){
  dir.create("Prelim_Figures")
}
figures <- "Prelim_Figures"

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
jpeg(filename = file.path(figures, "pH_SoilGrids.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(pH, par.settings = RdBuTheme, margin=FALSE)
dev.off()

OrgC <- raster(file.path(soil, "ORCDRC_weighted.tif"))
jpeg(filename = file.path(figures, "OrgC_SoilGrids.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(OrgC, par.settings = RdBuTheme, margin=FALSE)
dev.off()

clay <- raster(file.path(soil, "CLYPPT_weighted.tif"))
jpeg(filename = file.path(figures, "Clay_SoilGrids.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(clay, par.settings = RdBuTheme, margin=FALSE)
dev.off()

silt <- raster(file.path(soil, "SLTPPT_weighted.tif")) 
jpeg(filename = file.path(figures, "Silt_SoilGrids.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(silt, par.settings = RdBuTheme, margin=FALSE)
dev.off()

CEC <- raster(file.path(soil, "CECSOL_weighted.tif"))
jpeg(filename = file.path(figures, "CEC_SoilGrids.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(CEC, par.settings = RdBuTheme, margin=FALSE)
dev.off()

## Climate
chelsa <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\CHELSAData\\BioClim"

AnnPrecip <- raster(file.path(chelsa, "CHELSA_bio10_1.tif"))
jpeg(filename = file.path(figures, "AnnPrecip_Chelsa.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(AnnPrecip, par.settings = RdBuTheme, margin=FALSE)
dev.off()



PrecipSeason <- raster(file.path(chelsa, "CHELSA_bio10_15.tif"))
jpeg(filename = file.path(figures, "PrecipSeason_Chelsa.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(PrecipSeason, par.settings = RdBuTheme, margin=FALSE)
dev.off()



snow <- raster("I:\\sWorm\\Database\\snow\\snow_2015_sum.tif")
jpeg(filename = file.path(figures, "Snow.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(snow, par.settings = RdBuTheme, margin=FALSE)
dev.off()


PET <- raster("I:\\sWorm\\Database\\PET_he_annual\\pet_he_yr_TIF.tif")
jpeg(filename = file.path(figures, "PET.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
plot(PET)
dev.off()


PETSD <- raster("I:\\sWorm\\Database\\PET_he_monthly\\pet_he_SD.tif")
jpeg(filename = file.path(figures, "PETSD.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(PETSD, par.settings = RdBuTheme, margin=FALSE)
dev.off()

AI <- raster("I:\\sWorm\\Database\\AI_annual\\ai_yr_TIF.tif")
AI <- AI/10000 
jpeg(filename = file.path(figures, "AI.jpg"), width = 1500, height = 1000, quality = 150, pointsize = 14)
par(mar=c(5, 4, 1, 1))
levelplot(AI, par.settings = RdBuTheme, margin=FALSE)
dev.off()
