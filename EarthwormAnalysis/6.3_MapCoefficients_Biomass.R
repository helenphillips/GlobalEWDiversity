
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
  reg <- args[4] ## Which continent
  
  print(GLs_folder)
  print(models) 
  print(savefolder)
  print(reg)
  
  rasterOptions(tmpdir = "/work/phillips", chunksize = 524288, maxmemory = 134217728)
  
}

########################################################
# 2.5 Regions for analysis
########################################################

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


  if(!dir.exists(file.path(savefolder, reg))){
    dir.create(file.path(savefolder, reg))
  }
  #  data_out <- file.path(savefolder, reg)
  
  #################################################
  # 5. ABUNDANCE
  #################################################
  print("Creating biomass raster")
  print("Loading all rasters")
  
  print(file.path(GLs_folder, reg))
  
  bio12 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_12_BiomassCutScaled.tif"))
  bio15 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_15_BiomassCutScaled.tif"))
  
  snow <- raster(file.path(GLs_folder, reg, "Snow_newValues.tif"))

  ph <- raster(file.path(GLs_folder, reg,"PHIHOX_BiomassCutScaled.tif"))
  clay <- raster(file.path(GLs_folder, reg, "CLYPPT_BiomassCutScaled.tif"))
  silt <- raster(file.path(GLs_folder, reg, "SLTPPT_BiomassCutScaled.tif"))
  cec <- raster(file.path(GLs_folder, reg, "CECSOL_BiomassCutScaled.tif"))
  orgC <- raster(file.path(GLs_folder, reg, "ORCDRC_BiomassCutScaled.tif"))
  
  esa <- raster(file.path(GLs_folder, reg, "ESA_newValuesCropped.tif"))
  
  zeroLayer <- bio12
  zeroLayer[!is.na(zeroLayer)] <- 0
  
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
  snow4[snow4 == 4] <- round(fixedeffs['SnowMonths_cat4plus'], digits = 2)
  
    print("Adding together all Snow coefs....")
  AllSnow_coefs <- cover(zeroLayer, snow1, snow2, snow3, snow4,
                               filename = file.path(savefolder, reg, "biomass_allSnowCoefs.tif"), overwrite = TRUE)
  rm(snow1, snow2, snow3, snow4)
  
  ################## 
  print("Calculating interacting coefficients") 
  print("Soil....")
  
  phClay <- overlay(ph, clay, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "phClaybiomass.tif"), overwrite = TRUE)
  phClay <- calc(phClay, fun = function(x){round(x * fixedeffs['scalePH:scaleCLYPPT'], digits = 2)},
                 filename = file.path(savefolder, reg, "phClaybiomasscoef.tif"), overwrite = TRUE)
  
  phSilt <- overlay(ph, silt, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "phSiltbiomass.tif"), overwrite = TRUE)
  phSilt <- calc(phSilt, fun = function(x){round(x * fixedeffs['scalePH:scaleSLTPPT'], digits = 2)},
                 filename = file.path(savefolder, reg, "phSiltbiomasscoef.tif"), overwrite = TRUE)
  
  print("2 down...")
  
  phorgC <- overlay(ph, orgC, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "phorgCbiomass.tif"), overwrite = TRUE)
  phorgC <- calc(phorgC, fun = function(x){round(x * fixedeffs['scalePH:scaleORCDRC'], digits = 2)},
                 filename = file.path(savefolder, reg, "phorgCbiomasscoef.tif"), overwrite = TRUE)
  
  phCec <- overlay(ph, cec, fun = createInteractionCoef,
                   filename = file.path(savefolder, reg, "phCecbiomass.tif"), overwrite = TRUE)
  phCec <- calc(phCec, fun = function(x){round(x * fixedeffs['scalePH:scaleCECSOL'], digits = 2)},
                filename = file.path(savefolder, reg, "phCecbiomasscoef.tif"), overwrite = TRUE)
  
  print("4 down...")
  clayOrgC <- overlay(clay, orgC, fun = createInteractionCoef,
                   filename = file.path(savefolder, reg, "clayOrgCbiomass.tif"), overwrite = TRUE)
  clayOrgC <- calc(clayOrgC, fun = function(x){round(x * fixedeffs['scaleCLYPPT:scaleORCDRC'], digits = 2)},
                filename = file.path(savefolder, reg, "clayOrgcbiomasscoef.tif"), overwrite = TRUE)
  
  clayCec <- overlay(clay, cec, fun = createInteractionCoef,
                      filename = file.path(savefolder, reg, "clayCecbiomass.tif"), overwrite = TRUE)
  clayCec <- calc(clayCec, fun = function(x){round(x * fixedeffs['scaleCLYPPT:scaleCECSOL'], digits = 2)},
                   filename = file.path(savefolder, reg, "clayCecbiomasscoef.tif"), overwrite = TRUE)
  
  print("6 down, 2 to go....")
  siltCec <- overlay(silt, cec, fun = createInteractionCoef,
                     filename = file.path(savefolder, reg, "siltCecbiomass.tif"), overwrite = TRUE)
  siltCec <- calc(siltCec, fun = function(x){round(x * fixedeffs['scaleSLTPPT:scaleCECSOL'], digits = 2)},
                  filename = file.path(savefolder, reg, "siltCecbiomasscoef.tif"), overwrite = TRUE)
  
  orgCCec <- overlay(orgC, cec, fun = createInteractionCoef,
                     filename = file.path(savefolder, reg, "orgCCecbiomass.tif"), overwrite = TRUE)
  orgCCec <- calc(orgCCec, fun = function(x){round(x * fixedeffs['scaleORCDRC:scaleCECSOL'], digits = 2)},
                  filename = file.path(savefolder, reg, "orgCCecbiomasscoef.tif"), overwrite = TRUE)
  
  f_together <- function(a, b, c, d, e, f, g, h){
    round(a + b +c +d + e + f +g +h, digits = 2)
  }
  
  print("Adding together all Soil coefs....")
  AllSoil_coefs <- overlay(phClay,phSilt, phorgC, phCec, clayOrgC,siltCec, clayCec, orgCCec,
                               fun = f_together, filename = file.path(savefolder, reg, "biomass_allSoilCoefs.tif"), overwrite = TRUE)
  
  rm(phClay,phSilt, phorgC, phCec, clayOrgC,siltCec, clayCec, orgCCec)
  
  
  ##### Water retention
  print("Water retention....")

  
  
  clayBio12 <- overlay(clay, bio12, fun = createInteractionCoef,
                     filename = file.path(savefolder, reg, "clayBio12biomass.tif"), overwrite = TRUE)
  clayBio12 <- calc(clayBio12, fun = function(x){round(x * fixedeffs['scaleCLYPPT:bio10_12_scaled'], digits = 2)},
                  filename = file.path(savefolder, reg, "clayBio12biomasscoef.tif"), overwrite = TRUE)
  
  siltBio12 <- overlay(silt, bio12, fun = createInteractionCoef,
                       filename = file.path(savefolder, reg, "siltBio12biomass.tif"), overwrite = TRUE)
  siltBio12 <- calc(siltBio12, fun = function(x){round(x * fixedeffs['scaleSLTPPT:bio10_12_scaled'], digits = 2)},
                    filename = file.path(savefolder, reg, "siltBio12biomasscoef.tif"), overwrite = TRUE)
  
  clayBio15 <- overlay(clay, bio15, fun = createInteractionCoef,
                       filename = file.path(savefolder, reg, "clayBio15biomass.tif"), overwrite = TRUE)
  clayBio15 <- calc(clayBio15, fun = function(x){round(x * fixedeffs['scaleCLYPPT:bio10_15_scaled'], digits = 2)},
                    filename = file.path(savefolder, reg, "clayBio15biomasscoef.tif"), overwrite = TRUE)
  
  f_together <- function(a, b, c){
    round(a + b +c , digits = 2)
  }
  
  print("Adding together all water retention coefs....")
  AllWaterRetention_coefs <- overlay(clayBio12,siltBio12, clayBio15,
                           fun = f_together, filename = file.path(savefolder, reg, "biomass_allWaterRetentionCoefs.tif"), overwrite = TRUE)
  rm(clayBio12,siltBio12, clayBio15)
  #### Intercept
  ## Snow raster layer can become the intercept layer
  print("Creating intercept raster....")
  snow[!(is.na(snow))] <- fixedeffs['(Intercept)']
  intercept <- snow
  rm(snow)
  intercept <- writeRaster(intercept, filename=file.path(savefolder, reg, "intercept_biomasscoefs.tif"), format="GTiff", overwrite=TRUE)
  
  #### Remove all raster
  rm(bio12, bio15, ph, clay, silt, cec, orgC)
  
  ## Need to add this to something
  ######## Add them all together
  
  f_together <- function(a, b, c, d, e){
    round(a + b +c +d +e , digits = 2)
  }
  
  print("adding together all raster layers....")
  biomass_finalraster <- overlay(intercept, esa, AllSoil_coefs,
                                   AllWaterRetention_coefs, AllSnow_coefs, 
                                   fun = f_together, 
                                   filename = file.path(savefolder, reg, "BiomassFinalRaster.tif"), overwrite = TRUE)
  
  print("Done!")
