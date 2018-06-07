
########################################################
# 1. Load Libraries and Data
########################################################


library(raster)

########################################################
# 2. Set Working Directory or Cluster Info
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
  
  GLs_folder <- "I:\\sWorm\\ProcessedGLs\\Same_resolution\\regions"
  models <- "Models"
}else{ ## i.e. cluster
  args <- commandArgs(trailingOnly = TRUE)
  
  GLs_folder <- args[1] # GLs_dir
  models <- args[2] # models_dir
  savefolder <- args[3] # output_dir
  
  
  print(GLs_folder)
  print(models) 
  print(savefolder)
  
  rasterOptions(tmpdir = "/work/phillips", chunksize = 524288, maxmemory = 134217728)
  
}


#################################################
# 3. Load in models
#################################################
print("Loading in the biodiversity models")
load(file.path(models, "biomassmodel_full.rds"))

#################################################
# 4. Functions
#################################################

createCoef <- function(x, coef, ...){
  t <- x * coef
  round(t, digits = 2)
}

createInteractionCoef <- function(x, y){
  round(x * y, digits = 2)
}

#################################################
# START LOOP
#################################################

for(reg in regions){
  
  
  if(!dir.exists(file.path(savefolder, reg))){
    dir.create(file.path(savefolder, reg))
  }
  #  data_out <- file.path(savefolder, reg)
  
  #################################################
  # 5. ABUNDANCE
  #################################################
  print("Creating biomass raster")
  print("Loading all rasters")
  
  bio12 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_12_BiomassCutScaled.tif"))
  bio15 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_15_BiomassCutScaled.tif"))
  
  snow <- raster(file.path(GLs_folder, reg, "Snow_newValues"))
  
  ph <- raster(file.path(GLs_folder, reg,"PHIHOX_BiomassCutScaled.tif"))
  clay <- raster(file.path(GLs_folder, reg, "CLYPPT_BiomassCutScaled.tif"))
  silt <- raster(file.path(GLs_folder, reg, "SLTPPT_BiomassCutScaled.tif"))
  cec <- raster(file.path(GLs_folder, reg, "CECSOL_BiomassCutScaled.tif"))
  orgC <- raster(file.path(GLs_folder, reg, "ORCDRC_BiomassCutScaled.tif"))
  
  esa <- raster(file.path(GLs_folder, reg, "ESA_newValuesCropped.tif"))
  
  
  fixedeffs <- summary(biomass_model)$coefficients[,1]
  
  
  ##### ESA
  print("Habitat cover...")
  esa[esa == 0] <- NA ## They have an NA pixel number
  esa[esa == 60] <- 0 ## intercept is added later # Broadleaf deciduous
  print("2 done")
  esa[esa == 50] <- round(fixedeffs['ESABroadleaf evergreen forest'], digits = 2)
  esa[esa == 70] <- round(fixedeffs['ESANeedleleaf evergreen forest'], digits = 2)
  print("4 done")
  esa[esa == 90] <- round(fixedeffs['ESAMixed forest'], digits = 2)
  esa[esa == 110] <- round(fixedeffs['ESAHerbaceous with spare tree/shrub'], digits = 2)
  print("6 done")
  esa[esa == 120] <- round(fixedeffs['ESAShrub'], digits = 2)
  esa[esa == 130] <- round(fixedeffs['ESAHerbaceous'], digits = 2)
  print("8 done...last 3")
  esa[esa == 10] <- round(fixedeffs['ESAProduction - Herbaceous'], digits = 2)
  esa[esa == 12] <- round(fixedeffs['ESAProduction - Plantation'], digits = 2)
  esa[esa == 202] <- round(fixedeffs['ESABare area (unconsolidated'], digits = 2)
  
  ## Put all habitat covers where we don't have the habitat cover in the model to NA
  esa[esa > 10] <- NA
  
  print("Saving ESA layer")
  esa <- writeRaster(esa,  filename=file.path(savefolder, reg, "ESA_biomasscoefs.tif"), format="GTiff", overwrite=TRUE)
  
  ########### SNOW COVER
  print("Creating snow masks....")
  snow1 <- snow2 <- snow3 <-snow4 <- snow
  
  snow1[snow1 < 1 | snow1 > 1] <- NA
  snow1[snow1 == 1] <- round(fixedeffs['SnowMonths_cat1'], digits = 2)
  
  snow2[snow2 < 2 | snow2 > 2] <- NA
  snow2[snow2 == 2] <- round(fixedeffs['SnowMonths_cat2'], digits = 2)
  
  snow3[snow3 < 3 | snow3 > 3] <- NA
  snow3[snow3 == 3] <- round(fixedeffs['SnowMonths_cat3'], digits = 2)
  
  snow4[snow4 < 4 | snow4 > 4] <- NA
  snow4[snow4 == 4] <- round(fixedeffs['SnowMonths_cat4'], digits = 2)
  
  f_together <- function(a, b, c, d){
    round(a + b +c +d, digits = 2)
  }
  
  print("Adding together all Snow coefs....")
  AllBio1Snow_coefs <- overlay(snow1, snow2, snow3, snow4, fun = f_together, 
                               filename = file.path(savefolder, reg, "biomass_allSnowCoefs.tif"))
  rm(snow1, snow2, snow3, snow4)
  
  ################## 
    scaleCLYPPT:scaleORCDRC + scaleCLYPPT:scaleCECSOL + scaleSLTPPT:scaleCECSOL +  
    scaleORCDRC:scaleCECSOL 
  
  print("Calculating interacting coefficients") 
  print("Soil....")
  
  phClay <- overlay(ph, clay, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "phClaybiomass.tif"))
  phClay <- calc(phClay, fun = function(x){round(x * fixedeffs['scalePH:scaleCLYPPT'], digits = 2)},
                 filename = file.path(savefolder, reg, "phClaybiomass.tif"), overwrite = TRUE)
  
  phSilt <- overlay(ph, silt, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "phSiltbiomass.tif"))
  phSilt <- calc(phSilt, fun = function(x){round(x * fixedeffs['scalePH:scaleSLTPPT'], digits = 2)},
                 filename = file.path(savefolder, reg, "phSiltbiomass.tif"), overwrite = TRUE)
  
  phorgC <- overlay(ph, orgC, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "phorgCbiomass.tif"))
  phorgC <- calc(phorgC, fun = function(x){round(x * fixedeffs['scalePH:scaleORCDRC'], digits = 2)},
                 filename = file.path(savefolder, reg, "phorgCbiomass.tif"), overwrite = TRUE)
  
  phCec <- overlay(ph, cec, fun = createInteractionCoef,
                   filename = file.path(savefolder, reg, "phCecbiomass.tif"))
  phCec <- calc(phCec, fun = function(x){round(x * fixedeffs['scalePH:scaleCECSOL'], digits = 2)},
                filename = file.path(savefolder, reg, "phCecbiomass.tif"), overwrite = TRUE)
  
  
  + scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +  
    scaleCLYPPT:bio10_15_scaled