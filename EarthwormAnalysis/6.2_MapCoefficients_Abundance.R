
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

########################################################
# 2.5 Regions for analysis
########################################################
regions <- c("africa","asia","europe","latin_america","north_america","west_asia")

#################################################
# 3. Load in models
#################################################
print("Loading in the biodiversity models")
load(file.path(models, "abundancemodel_full.rds"))

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
  print("Creating abundance raster")
  print("Loading all rasters")
  
  
  bio1 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_1_AbundanceCutScaled.tif"))
  bio15 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_15_AbundanceCutScaled.tif"))
  
  snow <- raster(file.path(GLs_folder, reg, "Snow_newValues"))
  aridity <- raster(file.path(GLs_folder, reg, "Aridity_AbundanceScaled"))
  petsd <- raster(file.path(GLs_folder, reg, "PETSD_AbundanceScaled"))
  
  ph <- raster(file.path(GLs_folder, reg,"PHIHOX_AbundanceCutScaled.tif"))
  clay <- raster(file.path(GLs_folder, reg, "CLYPPT_AbundanceCutScaled.tif"))
  silt <- raster(file.path(GLs_folder, reg, "SLTPPT_AbundanceCutScaled.tif"))
  cec <- raster(file.path(GLs_folder, reg, "CECSOL_AbundanceCutScaled.tif"))
  orgC <- raster(file.path(GLs_folder, reg, "ORCDRC_AbundanceCutScaled.tif"))
  
  esa <- raster(file.path(GLs_folder, reg, "ESA_newValuesCropped.tif"))
  
  
  fixedeffs <- summary(abundance_model)$coefficients[,1]
  
  
  
  # fun <- function(x) { round(x * fixedeffs['scaleAridity'], digits = 2 )}
  # aridity <- calc(aridity, fun, filename = file.path(savefolder, "aridity_abundancecoef.tif"))
  # 
  
  # Intercept is added later, so can ignore 0
  
  # snow[snow == 1] <- round(fixedeffs['SnowMonths_cat1'], digits = 2)
  # snow[snow == 2] <- round(fixedeffs['SnowMonths_cat2'], digits = 2)
  # snow[snow == 3] <- round(fixedeffs['SnowMonths_cat3'], digits = 2)
  # snow[snow == 4] <- round(fixedeffs['SnowMonths_cat4plus'], digits = 2)
  # snow <- writeRaster(snow,  filename=file.path(savefolder, "snow_abundancecoefs.tif"), format="GTiff", overwrite=TRUE)
  
  # fun <- function(x) { round(x * fixedeffs['scalePH'], digits = 2 ) }
  # ph <- calc(ph, fun, filename = file.path(savefolder, "ph_abundancecoef.tif"))
  
  # fun <- function(x) { round(x * fixedeffs['scaleCLYPPT'], digits = 2 ) }
  # clay <- calc(clay, fun, filename = file.path(GLs_folder, "clay_abundancecoef.tif"))
  
  # fun <- function(x) { round(x * fixedeffs['scaleSLTPPT'], digits = 2 ) }
  # silt <- calc(silt, fun, filename = file.path(GLs_folder, "silt_abundancecoef.tif"))
  
  # fun <- function(x) { round(x * fixedeffs['scaleCECSOL'], digits = 2 ) }
  # cec <- calc(cec, fun, filename = file.path(savefolder, "cec_abundancecoef.tif"))
  
  # fun <- function(x) { round(x *fixedeffs['scaleORCDRC'], digits = 2 ) }
  # orgC <- calc(orgC, fun, filename = file.path(savefolder, "orgc_abundancecoef.tif"))
  
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
  esa[esa == 30] <- round(fixedeffs['ESACropland/Other vegetation mosaic'], digits = 2)
  
  ## Put all habitat covers where we don't have the habitat cover in the model to NA
  esa[esa > 10] <- NA
  
  print("Saving ESA layer")
  esa <- writeRaster(esa,  filename=file.path(savefolder, reg, "ESA_abundancecoefs.tif"), format="GTiff", overwrite=TRUE)
  #esa <- raster(file.path(GLs_folder, "ESA_coefs.tif"))
  #esa
  #str(esa)
  #summary(esa)
  #unique(getValues(esa))
  
  print("Calculating interacting coefficients") 
  print("Soil....")
  
  phClay <- overlay(ph, clay, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "phClayabundance.tif"))
  phClay <- calc(phClay, fun = function(x){round(x * fixedeffs['scalePH:scaleCLYPPT'], digits = 2)},
                 filename = file.path(savefolder, reg, "phClayabundancecoef.tif"), overwrite = TRUE)
  
  phCec <- overlay(ph, cec, fun = createInteractionCoef,
                   filename = file.path(savefolder, reg, "phCecabundance.tif"))
  phCec <- calc(phCec, fun = function(x){round(x * fixedeffs['scalePH:scaleCECSOL'], digits = 2)},
                filename = file.path(savefolder, reg, "phCecabundancecoef.tif"), overwrite = TRUE)
  
  clayCec <- overlay(clay, cec, fun = createInteractionCoef,
                     filename = file.path(savefolder, reg, "clayCecabundance.tif"))
  clayCec <- calc(clayCec, fun = function(x){round(x * fixedeffs['scaleCLYPPT:scaleCECSOL'], digits = 2)},
                  filename = file.path(savefolder, reg, "clayCecabundancecoef.tif"), overwrite = TRUE)
  
  siltCec <- overlay(silt, cec, fun = createInteractionCoef,
                     filename = file.path(savefolder, reg, "siltCecabundance.tif"))
  siltCec <- calc(siltCec, fun = function(x){round(x * fixedeffs['scaleSLTPPT:scaleCECSOL'], digits = 2)},
                  filename = file.path(savefolder, reg, "siltCecabundancecoef.tif"), overwrite = TRUE)
  
  cecOrgC <- overlay(cec, orgC, fun = createInteractionCoef,
                     filename = file.path(savefolder, reg, "cecCabundance.tif"))
  ceOrgC <- calc(cecOrgC, fun = function(x){round(x * fixedeffs['scaleCECSOL:scaleORCDRC'], digits = 2)},
                 filename = file.path(savefolder, reg, "cecCabundancecoef.tif"), overwrite = TRUE)
  
  
  
  
  f_together <- function(a, b, c, d, e){
    round(a + b +c +d +e, digits = 2)
  }
  
  print("Adding together all soil coefs....")
  allsoil_coefs <- overlay(phClay, phCec, clayCec, siltCec, cecOrgC, fun = f_together, 
                           filename = file.path(savefolder,  reg,"abundance_allsoilCoefs.tif"))
  
  rm(phClay)
  rm(phCec)
  rm(clayCec)
  rm(siltCec)
  rm(cecOrgC)
  
  print("Climate....")
  ## Creating snowmonth masks and stuff
  print("Creating snow masks....")
  snow0 <- snow1 <- snow2 <- snow3 <-snow4 <- snow
  snow0[snow0 > 0] <- NA
  snow0[snow0 == 0] <- 1 ## This just makes it easier for my multiplications
  
  snow1[snow1 < 1 | snow1 > 1] <- NA
  snow1mask <- snow1
  snow1[snow1 == 1] <- round(fixedeffs['SnowMonths_cat1'], digits = 2)
  
  snow2[snow2 < 2 | snow2 > 2] <- NA
  snow2mask <- snow2
  snow2mask[snow2mask == 2] <- 1
  snow2[snow2 == 2] <- round(fixedeffs['SnowMonths_cat2'], digits = 2)
  
  snow3[snow3 < 3 | snow3 > 3] <- NA
  snow3mask <- snow3
  snow3mask[snow3mask == 3] <- 1
  snow3[snow3 == 3] <- round(fixedeffs['SnowMonths_cat3'], digits = 2)
  
  snow4[snow4 < 4 | snow4 > 4] <- NA
  snow4mask <- snow4
  snow4mask[snow4mask == 4] <- 1
  snow4[snow4 == 4] <- round(fixedeffs['SnowMonths_cat4'], digits = 2)
  
  print("Interaction with snow....")
  ## Bio1 and Snow Months
  ## Snowmonth 0
  bio1snow0 <- overlay(bio1, snow0, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio1snow0abundance.tif"))
  bio1snow0 <- calc(bio1snow0, fun = function(x){round(x * fixedeffs['bio10_1_scaled'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio1snow0abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth1
  bio1snow1 <- overlay(bio1, snow1mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio1snow1abundance.tif"))
  bio1snow1 <- calc(bio1snow1, 
                    fun = function(x){round(x * fixedeffs['bio10_1_scaled:SnowMonths_cat1'] + fixedeffs['SnowMonths_cat1'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio1snow1abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth2
  bio1snow2 <- overlay(bio1, snow2mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio1snow2abundance.tif"))
  bio1snow2 <- calc(bio1snow2, 
                    fun = function(x){round(x * fixedeffs['bio10_1_scaled:SnowMonths_cat2'] + fixedeffs['SnowMonths_cat2'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio1snow2abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth3
  bio1snow3 <- overlay(bio1, snow3mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio1snow3abundance.tif"))
  bio1snow3 <- calc(bio1snow3, 
                    fun = function(x){round(x * fixedeffs['bio10_1_scaled:SnowMonths_cat3'] + fixedeffs['SnowMonths_cat3'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio1snow3abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth4
  bio1snow4 <- overlay(bio1, snow4mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio1snow4abundance.tif"))
  bio1snow4 <- calc(bio1snow4, 
                    fun = function(x){round(x * fixedeffs['bio10_1_scaled:SnowMonths_cat4'] + fixedeffs['SnowMonths_cat4'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio1snow4abundancecoef.tif"), overwrite = TRUE) 
  
  
  f_together <- function(a, b, c, d){
    round(a + b +c +d, digits = 2)
  }
  
  print("Adding together all bio1Snow coefs....")
  AllBio1Snow_coefs <- overlay(bio1snow1, bio1snow2, bio1snow3, bio1snow4, fun = f_together, 
                               filename = file.path(savefolder, reg, "abundance_allBio1SnowCoefs.tif"))
  rm(bio1snow1, bio1snow2, bio1snow3, bio1snow4)
  
  ## Bio15 and Snow Months
  ## Snowmonth 0
  bio15snow0 <- overlay(bio15, snow0, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "bio15snow0abundance.tif"))
  bio15snow0 <- calc(bio15snow0, fun = function(x){round(x * fixedeffs['bio10_15_scaled'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "bio15snow0abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth1
  bio15snow1 <- overlay(bio15, snow1mask, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "bio15snow1abundance.tif"))
  bio15snow1 <- calc(bio15snow1, 
                     fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat1'] + fixedeffs['SnowMonths_cat1'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "bio15snow1abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth2
  bio15snow2 <- overlay(bio15, snow2mask, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "bio15snow2abundance.tif"))
  bio15snow2 <- calc(bio15snow2, 
                     fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat2'] + fixedeffs['SnowMonths_cat2'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "bio15snow2abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth3
  bio15snow3 <- overlay(bio15, snow3mask, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "bio15snow3abundance.tif"))
  bio15snow3 <- calc(bio15snow3, 
                     fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat3'] + fixedeffs['SnowMonths_cat3'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "bio15snow3abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth4
  bio15snow4 <- overlay(bio15, snow4mask, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "bio15snow4abundance.tif"))
  bio15snow4 <- calc(bio15snow4, 
                     fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat4'] + fixedeffs['SnowMonths_cat4'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "bio15snow4abundancecoef.tif"), overwrite = TRUE) 
  
  print("Adding together all bio15Snow coefs....")
  AllBio15Snow_coefs <- overlay(bio15snow1, bio15snow2, bio15snow3, bio15snow4, fun = f_together, 
                                filename = file.path(savefolder, reg, "abundance_allBio15SnowCoefs.tif"))
  rm(bio15snow1, bio15snow2, bio15snow3, bio15snow4)
  
  ## SCalePETSD with snowmonths
  ## Snowmonth 0
  petsdsnow0 <- overlay(petsd, snow0, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "petsdsnow0abundance.tif"))
  petsdsnow0 <- calc(petsdsnow0, fun = function(x){round(x * fixedeffs['ScalePETSD'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "snow0petsd_abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth1
  petsdsnow1 <- overlay(petsd, snow1mask, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "petsdsnow1abundance.tif"))
  petsdsnow1 <- calc(petsdsnow1, 
                     fun = function(x){round(x * fixedeffs['ScalePETSD:SnowMonths_cat1'] + fixedeffs['SnowMonths_cat1'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "snow1petsd_abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth2
  petsdsnow2 <- overlay(petsd, snow2mask, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "petsdsnow2abundance.tif"))
  petsdsnow2 <- calc(petsdsnow2, 
                     fun = function(x){round(x * fixedeffs['ScalePETSD:SnowMonths_cat2'] + fixedeffs['SnowMonths_cat2'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "snow2petsd_abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth3
  petsdsnow3 <- overlay(petsd, snow3mask, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "petsdsnow3abundance.tif"))
  petsdsnow3 <- calc(petsdsnow3, 
                     fun = function(x){round(x * fixedeffs['ScalePETSD:SnowMonths_cat3'] + fixedeffs['SnowMonths_cat3'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "snow3petsd_abundancecoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth4
  petsdsnow4 <- overlay(petsd, snow4mask, fun = createInteractionCoef, 
                        filename = file.path(savefolder, reg, "petsdsnow4abundance.tif"))
  petsdsnow4 <- calc(petsdsnow4, 
                     fun = function(x){round(x * fixedeffs['ScalePETSD:SnowMonths_cat4'] + fixedeffs['SnowMonths_cat4'], digits = 2)}, 
                     filename = file.path(savefolder, reg, "snow4petsd_abundancecoef.tif"), overwrite = TRUE) 
  
  print("Adding together all petSDSnow coefs....")
  AllpetsdSnow_coefs <- overlay(petsdsnow1, petsdsnow2, petsdsnow3, petsdsnow4, fun = f_together, 
                                filename = file.path(savefolder, reg, "abundance_allpetSDSnowCoefs.tif"))
  rm(petsdsnow1, petsdsnow2, petsdsnow3, petsdsnow4)
  
  
  
  bio1bio15 <- overlay(bio1, bio15, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio1bio15abundance.tif"))
  bio1bio15 <- calc(bio1bio15, fun = function(x){round(x * fixedeffs['bio10_1_scaled:bio10_15_scaled'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio1bio15abundancecoef.tif"), overwrite = TRUE) 
  
  bio1arid <- overlay(bio1, aridity, fun = createInteractionCoef, 
                      filename = file.path(savefolder, reg, "bio1bio15abundance.tif"))
  bio1arid <- calc(bio1arid, fun = function(x){round(x * fixedeffs['bio10_1_scaled:scaleAridity'], digits = 2)}, 
                   filename = file.path(savefolder, reg, "bio1aridabundancecoef.tif"), overwrite = TRUE) 
  
  bio1petsd <- overlay(bio1, petsd, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio1petsdabundance.tif"))
  bio1petsd <- calc(bio1petsd, fun = function(x){round(x * fixedeffs['bio10_1_scaled:ScalePETSD'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio1petsdabundancecoef.tif"), overwrite = TRUE) 
  
  
  
  f_together <- function(a, b, c){
    round(a + b +c, digits = 2)
  }
  
  print("Adding together all other climate coefs....")
  allclimate_coefs <- overlay(bio1bio15, bio1arid, bio1petsd, fun = f_together, 
                              filename = file.path(savefolder, reg, "abundance_allotherClimateCoefs.tif"))
  
  
  
  
  print("Water retention....")
  claybio15 <- overlay(bio15, clay, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "claybio15abundance.tif"))
  claybio15 <- calc(claybio15, fun = function(x){round(x * fixedeffs['scaleCLYPPT:bio10_15_scaled'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "claybio15abundancecoef.tif"), overwrite = TRUE) 
  
  claypetsd <- overlay(clay, petsd, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "petsdclayabundance.tif"))
  claypetsd <- calc(claypetsd, fun = function(x){round(x * fixedeffs['scaleCLYPPT:ScalePETSD'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "petsdclayabundancecoef.tif"), overwrite = TRUE) 
  
  f_together <- function(a, b){
    round(a + b, digits = 2)
  }
  
  print("Adding together all water retention coefs....")
  allwaterretention_coefs <- overlay(claypetsd, claybio15, fun = f_together, 
                                     filename = file.path(savefolder, reg, "abundance_allWaterRetentionCoefs.tif"))
  
  
  #### Intercept
  ## Snow raster layer can become the intercept layer
  print("Creating intercept raster....")
  snow[!(is.na(snow))] <- fixedeffs['(Intercept)']
  intercept <- snow
  rm(snow)
  intercept <- writeRaster(intercept, filename=file.path(savefolder, reg, "intercept_abundancecoefs.tif"), format="GTiff", overwrite=TRUE)
  
  #### Remove all raster
  rm(bio1, bio15, aridity, petsd, ph, clay, silt, cec, org)
  
  ## Need to add this to something
  ######## Add them all together
  
  f_together <- function(a, b, c, d, e, f, g, h){
    round(a + b +c +d +e +f +g + h, digits = 2)
  }
  
  print("adding together all raster layers....")
  abundance_finalraster <- overlay(intercept, esa, allsoil_coefs, AllBio1Snow_coefs, AllBio15Snow_coefs,
                                   AllpetsdSnow_coefs
                                   allwaterretention_coefs, allclimate_coefs, 
                                   fun = f_together, 
                                   filename = file.path(savefolder, reg, "AbundanceFinalRaster.tif"))
  
  print("Done!")
}