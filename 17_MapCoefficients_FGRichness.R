
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


#################################################
# 3. Load in models
#################################################
print("Loading in the biodiversity models")
load(file.path(models, "fgrichnessmodel.rds"))

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

  if(!dir.exists(file.path(savefolder, reg))){
    dir.create(file.path(savefolder, reg))
  }
 #  data_out <- file.path(savefolder, reg)

  #################################################
  # 5. FG Richness
  #################################################
  print("Creating Functional group richness raster")
  print("Loading all rasters")
  
  print(file.path(GLs_folder, reg))
  
  print("Chelsa data")
  bio4 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_4_FGRichnessCutScaled.tif"))
  bio15 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_15_FGRichnessCutScaled.tif"))
  
  print("Other climate data")
  snow <- raster(file.path(GLs_folder, reg, "Snow_newValues.tif"))
  aridity <- raster(file.path(GLs_folder, reg, "Aridity_FGRichnessScaled.tif"))
  pet <- raster(file.path(GLs_folder, reg, "PETyr_FGRichnessScaled.tif"))
  
  print("Soil data")
  silt <- raster(file.path(GLs_folder, reg, "SLTPPT_FGRichnessCutScaled.tif"))
  orgC <- raster(file.path(GLs_folder, reg, "ORCDRC_FGRichnessCutScaled.tif"))

  fixedeffs <- summary(fgrichness_model)$coefficients[,1]
  
  #############################################################
  
  print("Calculating interacting coefficients") 
  print("Soil....")
  
  
  siltOrgC <- overlay(silt, orgC, fun = createInteractionCoef,
                     filename = file.path(savefolder, reg, "siltOrgCFGRichness.tif"), overwrite = TRUE)
  siltOrgC <- calc(siltOrgC, fun = function(x){round(x * fixedeffs['scaleSLTPPT:scaleORCDRC'], digits = 2)},
                  filename = file.path(savefolder, reg, "siltOrgCFGRichnesscoef.tif"), overwrite = TRUE)
  
 
  
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
  snow4[snow4 == 4] <- round(fixedeffs['SnowMonths_cat4plus'], digits = 2)
  
  print("Interaction with snow....")
  ## Bio1 and Snow Months
  ## Snowmonth 0
  bio15snow0 <- overlay(bio15, snow0, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow0fgrichness.tif"), overwrite = TRUE)
  bio15snow0 <- calc(bio15snow0, fun = function(x){round(x * fixedeffs['bio10_15_scaled'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15snow0fgrichnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth1
  bio15snow1 <- overlay(bio15, snow1mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow1fgrichness.tif"), overwrite = TRUE)
  bio15snow1 <- calc(bio15snow1, 
                    fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat1'] + fixedeffs['SnowMonths_cat1'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15snow1fgrichnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth2
  bio15snow2 <- overlay(bio15, snow2mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow2fgrichness.tif"), overwrite = TRUE)
  bio15snow2 <- calc(bio15snow2, 
                    fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat2'] + fixedeffs['SnowMonths_cat2'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15snow2fgrichnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth3
  bio15snow3 <- overlay(bio15, snow3mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow3fgrichness.tif"), overwrite = TRUE)
  bio15snow3 <- calc(bio15snow3, 
                    fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat3'] + fixedeffs['SnowMonths_cat3'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15snow3fgrichnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth4
  bio15snow4 <- overlay(bio15, snow4mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow4fgrichness.tif"), overwrite = TRUE)
  bio15snow4 <- calc(bio15snow4, 
                    fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat4plus'] + fixedeffs['SnowMonths_cat4plus'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio51snow4fgrichnesscoef.tif"), overwrite = TRUE) 
  

  print("Adding together all bio1Snow coefs....")
  AllBio15Snow_coefs <- cover(bio15snow0, bio15snow1, bio15snow2, bio15snow3, bio15snow4, 
                               filename = file.path(savefolder, reg, "fgrichness_allBio15SnowCoefs.tif"), overwrite = TRUE)
  rm(bio15snow0, bio15snow1, bio15snow2, bio15snow3, bio15snow4)
  
 
  print("Interacting climate variables....")
  
  bio4bio15 <- overlay(bio4, bio15, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio4bio15fgrichness.tif"), overwrite = TRUE)
  bio4bio15 <- calc(bio4bio15, fun = function(x){round(x * fixedeffs['bio10_4_scaled:bio10_15_scaled'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio4bio15abundancecoef.tif"), overwrite = TRUE) 
  
  bio15arid <- overlay(bio15, aridity, fun = createInteractionCoef, 
                      filename = file.path(savefolder, reg, "aridbio15abundance.tif"), overwrite = TRUE)
  bio15arid <- calc(bio1arid, fun = function(x){round(x * fixedeffs['bio10_15_scaled:scaleAridity'], digits = 2)}, 
                   filename = file.path(savefolder, reg, "bio15aridfgrichnesscoef.tif"), overwrite = TRUE) 
  
  bio15pet <- overlay(bio15, pet, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15petfgrichness.tif"), overwrite = TRUE)
  bio15pet <- calc(bio15pet, fun = function(x){round(x * fixedeffs['bio10_15_scaled:ScalePET'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15petfgrichnesscoef.tif"), overwrite = TRUE) 
  
  
  
  f_together <- function(a, b, c){
    round(a + b +c, digits = 2)
  }
  
  print("Adding together all other climate coefs....")
  allclimate_coefs <- overlay(bio4bio15, bio15arid, bio15pet, fun = f_together, 
                              filename = file.path(savefolder, reg, "fgrichness_allotherClimateCoefs.tif"), overwrite = TRUE)
  
  
  
  
  print("Water retention....")
  siltpet <- overlay(silt, pet, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "petsiltfgrichness.tif"), overwrite = TRUE)
  siltpet <- calc(siltpet, fun = function(x){round(x * fixedeffs['scaleSLTPPT:ScalePET'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "petsiltfgrichnesscoef.tif"), overwrite = TRUE) 
  
  #### Intercept
  ## Snow raster layer can become the intercept layer
  print("Creating intercept raster....")
  snow[!(is.na(snow))] <- fixedeffs['(Intercept)']
  intercept <- snow
  rm(snow)
  intercept <- writeRaster(intercept, filename=file.path(savefolder, reg, "intercept_fgrichnesscoefs.tif"), format="GTiff", overwrite=TRUE)
  
  #### Remove all raster
  rm(bio4, bio15, aridity, pet, silt, orgC)
  
  ## Need to add this to something
  ######## Add them all together
  
  f_together <- function(a, b, c, d, e){
    round(a + b +c +d +e , digits = 2)
  }
  
  print("adding together all raster layers....")
  fgrichness_finalraster <- overlay(intercept, 
                                   siltpet,
                                   siltOrgC,
                                   AllBio15Snow_coefs,
                                   allclimate_coefs,
                                   fun = f_together, 
                                   filename = file.path(savefolder, reg, "FGRichnessFinalRaster.tif"), overwrite = TRUE)
  
  print("Done!")
