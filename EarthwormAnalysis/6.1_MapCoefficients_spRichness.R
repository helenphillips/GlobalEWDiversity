
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
load(file.path(models, "richnessmodel_full.rds"))

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
# 5. RICHNESS
#################################################
  print("Creating richness raster")
  print("Loading all rasters")
  
  bio4 <- raster(file.path(GLs_folder,reg, "CHELSA_bio10_4_RichnessCutScaled.tif"))
  bio15 <- raster(file.path(GLs_folder,reg, "CHELSA_bio10_15_RichnessCutScaled.tif"))

  snow <- raster(file.path(GLs_folder,reg, "Snow_newValues"))
  aridity <- raster(file.path(GLs_folder,reg, "Aridity_RichnessScaled"))
  pet <- raster(file.path(GLs_folder,reg, "PETyr_RichnessScaled"))

  ph <- raster(file.path(GLs_folder,reg,"PHIHOX_RichnessCutScaled.tif"))
  clay <- raster(file.path(soil_GLs,reg,"CLYPPT_RichnessCutScaled.tif"))
  silt <- raster(file.path(soil_GLs,reg,"SLTPPT_RichnessCutScaled.tif"))
  cec <- raster(file.path(GLs_folder,reg,"CECSOL_RichnessCutScaled.tif"))
  orgC <- raster(file.path(GLs_folder,reg,"ORCDRC_RichnessCutScaled.tif"))

  esa <- raster(file.path(GLs_folder,reg,"ESA_newValuesCropped.tif"))


  fixedeffs <- summary(richness_model)$coefficients[,1]


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
  
  ## Put all habitat covers where we don't have the habitat cover in the model to NA
  esa[esa > 10] <- NA
  
  print("Saving ESA layer")
  esa <- writeRaster(esa,  filename=file.path(savefolder, "ESA_richnesscoefs.tif"), format="GTiff", overwrite=TRUE)
  #esa <- raster(file.path(GLs_folder, "ESA_coefs.tif"))

  print("Calculating interacting coefficients") 
  print("Soil....")

  
  phCec <- overlay(ph, cec, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "phCecrichness.tif"))
  phCec <- calc(phCec, fun = function(x){round(x * fixedeffs['scalePH:scaleCECSOL'], digits = 2)},
                 filename = file.path(savefolder, reg, "phCecrichnesscoef.tif"), overwrite = TRUE)
  
  phOrgC <- overlay(ph, orgC, fun = createInteractionCoef,
                   filename = file.path(savefolder, reg, "phOrgCrichness.tif"))
  phOrgC <- calc(phOrgC, fun = function(x){round(x * fixedeffs['scalePH:scaleORCDRC'], digits = 2)},
                filename = file.path(savefolder, reg, "phOrgCrichnesscoef.tif"), overwrite = TRUE)
  
  claycec <- overlay(clay, cec, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "clayCecrichness.tif"))
  claycec <- calc(claycec, fun = function(x){round(x * fixedeffs['scaleCLYPPT:scaleCECSOL'], digits = 2)},
                 filename = file.path(savefolder, reg, "clayCecrichnesscoef.tif"), overwrite = TRUE)
  
  siltOrgC <- overlay(silt, orgC, fun = createInteractionCoef,
                    filename = file.path(savefolder, reg, "siltOrgCrichness.tif"))
  siltOrgC <- calc(siltOrgC, fun = function(x){round(x * fixedeffs['scaleSLTPPT:scaleORCDRC'], digits = 2)},
                 filename = file.path(savefolder, reg, "siltOrgCrichnesscoef.tif"), overwrite = TRUE)
  
  cecOrgC <- overlay(cec, orgC, fun = createInteractionCoef,
                      filename = file.path(savefolder, reg, "cecOrgCrichness.tif"))
  cecOrgC <- calc(cecOrgC, fun = function(x){round(x * fixedeffs['scaleCECSOL:scaleORCDRC'], digits = 2)},
                   filename = file.path(savefolder, reg, "cecOrgCrichnesscoef.tif"), overwrite = TRUE)
  
  f_together <- function(a, b, c, d, e){
    round(a + b +c +d + e, digits = 2)
  }
  
  print("Adding together all soil coefs....")
  Allsoil_coefs <- overlay(phCec,phOrgC, claycec,siltOrgC, cecOrgC, 
                          fun = f_together, filename = file.path(savefolder, reg, "richness_allSoilCoefs.tif"))
  rm(phCec,phOrgC, claycec,siltOrgC, cecOrgC)
  
  #################################
  print("Climate....")
  
  bio4bio15 <- overlay(bio4, bio15, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio4bio15richness.tif"))
  bio4bio15 <- calc(bio4bio15, fun = function(x){round(x * fixedeffs['bio10_4_scaled:bio10_15_scaled'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio4bio15richnesscoef.tif"), overwrite = TRUE) 
  
  
  aridbio15 <- overlay(aridity, bio15, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "aridbio15richness.tif"))
  aridbio15 <- calc(aridbio15, fun = function(x){round(x * fixedeffs['bio10_15_scaled:scaleAridity'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "aridbio15richnesscoef.tif"), overwrite = TRUE) 
  
  
  petbio15 <- overlay(pet, bio15, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "petbio15richness.tif"))
  petbio15 <- calc(petbio15, fun = function(x){round(x * fixedeffs['bio10_15_scaled:ScalePET'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "petbio15richnesscoef.tif"), overwrite = TRUE) 
  
  f_together <- function(a, b, c){
    round(a + b +c , digits = 2)
  }
  
  print("Adding together all climate coefs....")
  Allclimate_coefs <- overlay(bio4bio15,aridbio15, petbio15,
                           fun = f_together, filename = file.path(savefolder, reg, "richness_allclimateCoefs.tif"))
  
  rm(bio4bio15,aridbio15, petbio15)
  
  #######################
  ## Snow months
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
  bio4snow0 <- overlay(bio4, snow0, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio1snow0richness.tif"))
  bio4snow0 <- calc(bio4snow0, fun = function(x){round(x * fixedeffs['bio10_4_scaled'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio4snow0richnesscoef.tif"), overwrite = TRUE) 
  
  
  ## Snowmonth1
  bio4snow1 <- overlay(bio4, snow1mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio4snow1richness.tif"))
  bio4snow1 <- calc(bio4snow1, 
                    fun = function(x){round(x * fixedeffs['bio10_4_scaled:SnowMonths_cat1'] + fixedeffs['SnowMonths_cat1'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio4snow1richnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth2
  bio4snow2 <- overlay(bio4, snow2mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio4snow2richness.tif"))
  bio4snow2 <- calc(bio4snow2, 
                    fun = function(x){round(x * fixedeffs['bio10_4_scaled:SnowMonths_cat2'] + fixedeffs['SnowMonths_cat2'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio4snow2richnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth3
  bio4snow3 <- overlay(bio4, snow3mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio4snow3richness.tif"))
  bio4snow3 <- calc(bio4snow3, 
                    fun = function(x){round(x * fixedeffs['bio10_4_scaled:SnowMonths_cat3'] + fixedeffs['SnowMonths_cat3'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio4snow3richnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth4
  bio4snow4 <- overlay(bio4, snow4mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio4snow4richness.tif"))
  bio4snow4 <- calc(bio4snow4, 
                    fun = function(x){round(x * fixedeffs['bio10_4_scaled:SnowMonths_cat4plus'] + fixedeffs['SnowMonths_cat4plus'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio4snow4richnesscoef.tif"), overwrite = TRUE) 
  
  
  print("Adding together all bio4Snow coefs....")
  f_together <- function(a, b, c, d){
    round(a + b +c + d , digits = 2)
  }
  Allbio4Snow_coefs <- overlay(bio4snow1, bio4snow2, bio4snow3, bio4snow4, fun = f_together, 
                                filename = file.path(savefolder, reg, "richness_allbio4nowCoefs.tif"))
  rm(bio4snow1, bio4snow2, bio4snow3, bio4snow4)
  
  
  ## Bio15 and Snow Months
  ## Snowmonth 0
  bio15snow0 <- overlay(bio15, snow0, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow0richness.tif"))
  bio15snow0 <- calc(bio15snow0, fun = function(x){round(x * fixedeffs['bio10_15_scaled'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15snow0richnesscoef.tif"), overwrite = TRUE) 
  
  
  ## Snowmonth1
  bio15snow1 <- overlay(bio15, snow1mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow1richness.tif"))
  bio15snow1 <- calc(bio15snow1, 
                    fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat1'] + fixedeffs['SnowMonths_cat1'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15snow1richnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth2
  bio15snow2 <- overlay(bio15, snow2mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow2richness.tif"))
  bio15snow2 <- calc(bio15snow2, 
                    fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat2'] + fixedeffs['SnowMonths_cat2'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15snow2richnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth3
  bio15snow3 <- overlay(bio15, snow3mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow3richness.tif"))
  bio15snow3 <- calc(bio15snow3, 
                    fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat3'] + fixedeffs['SnowMonths_cat3'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15snow3richnesscoef.tif"), overwrite = TRUE) 
  
  ## Snowmonth4
  bio15snow4 <- overlay(bio15, snow4mask, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "bio15snow4richness.tif"))
  bio15snow4 <- calc(bio15snow4, 
                    fun = function(x){round(x * fixedeffs['bio10_15_scaled:SnowMonths_cat4plus'] + fixedeffs['SnowMonths_cat4plus'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "bio15snow4richnesscoef.tif"), overwrite = TRUE) 
  print("Adding together all bio15Snow coefs....")
  f_together <- function(a, b, c, d){
    round(a + b +c + d , digits = 2)
  }
  Allbio15Snow_coefs <- overlay(bio15snow1, bio15snow2, bio15snow3, bio15snow4, fun = f_together, 
                               filename = file.path(savefolder, reg, "richness_allbio15SnowCoefs.tif"))
  rm(bio4snow1, bio4snow2, bio4snow3, bio4snow4)
  
  
  ######################################
  #
  print("Water retention.....")
  
     + scaleSLTPPT:scaleAridity
  
  claybio15 <- overlay(clay, bio15, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "claybio15richness.tif"))
  claybio15 <- calc(claybio15, fun = function(x){round(x * fixedeffs['scaleCLYPPT:bio10_15_scaled'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "claybio15richnesscoef.tif"), overwrite = TRUE) 
  
  siltpet <- overlay(silt, pet, fun = createInteractionCoef, 
                       filename = file.path(savefolder, reg, "siltpetrichness.tif"))
  siltpet <- calc(siltpet, fun = function(x){round(x * fixedeffs['scaleSLTPPT:ScalePET'], digits = 2)}, 
                    filename = file.path(savefolder, reg, "siltpetrichnesscoef.tif"), overwrite = TRUE) 
  
  siltarid <- overlay(silt, aridity, fun = createInteractionCoef, 
                     filename = file.path(savefolder, reg, "siltaridrichness.tif"))
  siltarid <- calc(siltarid, fun = function(x){round(x * fixedeffs['scaleSLTPPT:scaleAridity'], digits = 2)}, 
                  filename = file.path(savefolder, reg, "siltaridrichnesscoef.tif"), overwrite = TRUE) 
  
  f_together <- function(a, b, c){
    round(a + b +c, digits = 2)
  }
  AllWaterRetention_coefs <- overlay(claybio15, siltpet, siltarid, fun = f_together, 
                                filename = file.path(savefolder, reg, "richness_allwaterretentionCoefs.tif"))
  
  rm(claybio15, siltpet, siltarid)
  
  
  ####################################
  #### Intercept
  ## Snow raster layer can become the intercept layer
  print("Creating intercept raster....")
  snow[!(is.na(snow))] <- fixedeffs['(Intercept)']
  intercept <- snow
  rm(snow)
  intercept <- writeRaster(intercept, filename=file.path(savefolder, reg, "intercept_richnesscoefs.tif"), format="GTiff", overwrite=TRUE)
  
  #### Remove all raster
  rm(bio4, bio15, aridity, pet, ph, clay, silt, cec, orgC)
  
  
  #########################
  print("Adding all rasters together")
  
  f_together <- function(a, b, c, d, e, f, g){
    round(a + b +c +d +e + f + g, digits = 2)
  }
  
  
  spR_finalraster <- overlay(esa, Allsoil_coefs, Allclimate_coefs, Allbio4Snow_coefs, Allbio15Snow_coefs,
                             AllWaterRetention_coefs, intercept,
                             fun = f_together, 
                             filename = file.path(savefolder, reg, "spRFinalRaster.tif"))
  print("Done!") 
}
