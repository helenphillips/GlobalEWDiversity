
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
  
  print("Habitat cover")
  esa <- raster(file.path(GLs_folder, reg, "ESA_newValuesCropped.tif"))
  
  print("Chelsa data")
  bio1 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_1_FGRichnessCutScaled.tif"))
  bio15 <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_15_FGRichnessCutScaled.tif"))
  
  print("Other climate data")
  # snow <- raster(file.path(GLs_folder, reg, "Snow_newValues_WGS84.tif"))
  aridity <- raster(file.path(GLs_folder, reg, "Aridity_FGRichnessScaled.tif"))
  petsd <- raster(file.path(GLs_folder, reg, "PETSD_FGRichnessScaled.tif"))
  
  print("Soil data")
  clay <- raster(file.path(GLs_folder, reg, "CLYPPT_FGRichnessCutScaled.tif"))
  ph <- raster(file.path(GLs_folder, reg, "PHIHOX_FGRichnessCutScaled.tif"))
  silt <- raster(file.path(GLs_folder, reg, "SLTPPT_FGRichnessCutScaled.tif"))

    fixedeffs <- summary(fgrichness_model)$coefficients[,1]
  
  #############################################################
  print("Doing the ESA layer...")
  
  esa[esa == 0] <- NA ## They have an NA pixel number 
  esa[esa == 60] <- 0 ## intercept is added later
  print("2 done")
  esa[esa == 50] <- round(fixedeffs['ESABroadleaf evergreen forest'], digits = 2)
  esa[esa == 70] <- round(fixedeffs['ESANeedleleaf evergreen forest'], digits = 2)
  print("4 done")
  esa[esa == 90] <- round(fixedeffs['ESAMixed forest'], digits = 2)
  esa[esa == 110] <- round(fixedeffs['ESAHerbaceous with spare tree/shrub'], digits = 2)
  print("6 done")
  esa[esa == 12] <- round(fixedeffs['ESAProduction - Plantation'], digits = 2)
  esa[esa == 120] <- round(fixedeffs['ESAShrub'], digits = 2)
  print("8 done...last 2")
  esa[esa == 130] <- round(fixedeffs['ESAHerbaceous'], digits = 2)
  esa[esa == 10] <- round(fixedeffs['ESAProduction - Herbaceous'], digits = 2)
  print("all done")
  
    ## Put all habitat covers where we don't have the habitat cover in the model to NA
  esa[esa > 10] <- NA
  
  print("Saving ESA layer")
  esa <- writeRaster(esa,  filename=file.path(savefolder, "ESA_FGRichnesscoefs.tif"), format="GTiff", overwrite=TRUE)


  print("Soil....no interactions. Only single effects")
  
  
  
  
 
  print("Climate....")
  print("Single effects....")
  
  print("Interacting climate variables....")
 
  bio1arid <- overlay(bio1, aridity, fun = createInteractionCoef, 
                      filename = file.path(savefolder, reg, "aridbio1fgrichnesscoef.tif"), overwrite = TRUE)
  bio1arid <- calc(bio1arid, fun = function(x){round(x * fixedeffs['bio10_1_scaled:scaleAridity'], digits = 2)}, 
                   filename = file.path(savefolder, reg, "bio1aridfgrichnesscoef.tif"), overwrite = TRUE) 
  
  bio1petsd <- overlay(bio1, petsd, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio1petsdfgrichness.tif"), overwrite = TRUE)
  bio1petsd <- calc(bio1petsd, fun = function(x){round(x * fixedeffs['bio10_1_scaled:ScalePETSD'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio1petsdfgrichnesscoef.tif"), overwrite = TRUE) 
  
  aridpetsd <- overlay(aridity, petsd, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "aridpetsdfgrichness.tif"), overwrite = TRUE)
  aridpetsd <- calc(aridpetsd, fun = function(x){round(x * fixedeffs['scaleAridity:ScalePETSD'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "aridpetsdfgrichnesscoef.tif"), overwrite = TRUE) 
  
  
  
  f_together <- function(a, b, c){
    round(a + b +c, digits = 2)
  }
  
  print("Adding together all other climate coefs....")
  allclimate_coefs <- overlay(bio1arid, bio1petsd, aridpetsd, fun = f_together, 
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
