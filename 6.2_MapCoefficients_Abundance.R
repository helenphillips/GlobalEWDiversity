
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
# 5. ABUNDANCE
#################################################
print("Creating abundance raster")
print("Loading all rasters")


bio1 <- raster(file.path(GLs_folder, "CHELSA_bio10_1_AbundanceCutScaled.tif"))
bio15 <- raster(file.path(GLs_folder, "CHELSA_bio10_15_AbundanceCutScaled.tif"))

snow <- raster(file.path(GLs_folder, "Snow_newValues"))
aridity <- raster(file.path(GLs_folder, "Aridity_AbundanceScaled"))
petsd <- raster(file.path(GLs_folder, "PETSD_AbundanceScaled"))

ph <- raster(file.path(GLs_folder,"PHIHOX_AbundanceCutScaled.tif"))
clay <- raster(file.path(soil_GLs,"CLYPPT_AbundanceCutScaled.tif"))
silt <- raster(file.path(soil_GLs,"SLTPPT_AbundanceCutScaled.tif"))
cec <- raster(file.path(GLs_folder,"CECSOL_AbundanceCutScaled.tif"))
orgC <- raster(file.path(GLs_folder,"ORCDRC_AbundanceCutScaled.tif"))

esa <- raster(file.path(GLs_folder,"ESA_newValuesCropped.tif"))


fixedeffs <- summary(richness_model)$coefficients[,1]

print("Calculating individual coefficients")
print("Climate coefficients....")

fun <- function(x) { round(x * fixedeffs['bio10_1_scaled'], digits = 2 )}
bio1 <- calc(bio7, fun, filename = file.path(savefolder, "bio10_1_abundancecoef.tif"))

fun <- function(x) { round(x * fixedeffs['bio10_15_scaled'], digits = 2 )}
bio15 <- calc(bio15, fun, filename = file.path(savefolder, "bio10_15_abundancecoef.tif"))

fun <- function(x) { round(x * fixedeffs['scaleAridity'], digits = 2 )}
aridity <- calc(aridity, fun, filename = file.path(savefolder, "aridity_abundancecoef.tif"))

fun <- function(x) { round(x * fixedeffs['ScalePETSD'], digits = 2 )}
petsd <- calc(petsd, fun, filename = file.path(savefolder, "petsd_abundancecoef.tif"))

# Intercept is added later, so can ignore 0

snow[snow == 1] <- round(fixedeffs['SnowMonths_cat1'], digits = 2)
snow[snow == 2] <- round(fixedeffs['SnowMonths_cat2'], digits = 2)
snow[snow == 3] <- round(fixedeffs['SnowMonths_cat3'], digits = 2)
snow[snow == 4] <- round(fixedeffs['SnowMonths_cat4plus'], digits = 2)
snow <- writeRaster(snow,  filename=file.path(savefolder, "snow_abundancecoefs.tif"), format="GTiff", overwrite=TRUE)

print("Soil coefficients....")
fun <- function(x) { round(x * fixedeffs['scalePH'], digits = 2 ) }
ph <- calc(ph, fun, filename = file.path(savefolder, "ph_abundancecoef.tif"))

fun <- function(x) { round(x * fixedeffs['scaleCLYPPT'], digits = 2 ) }
clay <- calc(clay, fun, filename = file.path(GLs_folder, "clay_abundancecoef.tif"))

fun <- function(x) { round(x * fixedeffs['scaleSLTPPT'], digits = 2 ) }
silt <- calc(silt, fun, filename = file.path(GLs_folder, "silt_abundancecoef.tif"))

fun <- function(x) { round(x * fixedeffs['scaleCECSOL'], digits = 2 ) }
cec <- calc(cec, fun, filename = file.path(savefolder, "cec_abundancecoef.tif"))

fun <- function(x) { round(x *fixedeffs['scaleORCDRC'], digits = 2 ) }
orgC <- calc(orgC, fun, filename = file.path(savefolder, "orgc_abundancecoef.tif"))

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
esa <- writeRaster(esa,  filename=file.path(savefolder, "ESA_abundancecoefs.tif"), format="GTiff", overwrite=TRUE)
#esa <- raster(file.path(GLs_folder, "ESA_coefs.tif"))
#esa
#str(esa)
#summary(esa)
#unique(getValues(esa))

print("Calculating interacting coefficients") 
print("Soil....")

phClay <- overlay(ph, clay, fun = createInteractionCoef,
                  filename = file.path(savefolder, "phClayabundance.tif"))
phClay <- calc(phClay, fun = function(x){round(x * fixedeffs['scalePH:scaleCLYPPT'], digits = 2)},
               filename = file.path(savefolder, "phClayabundance.tif"), overwrite = TRUE)

phCec <- overlay(ph, cec, fun = createInteractionCoef,
                  filename = file.path(savefolder, "phCecabundance.tif"))
phCec <- calc(phCec, fun = function(x){round(x * fixedeffs['scalePH:scaleCECSOL'], digits = 2)},
               filename = file.path(savefolder, "phCecabundance.tif"), overwrite = TRUE)

clayCec <- overlay(clay, cec, fun = createInteractionCoef,
                 filename = file.path(savefolder, "clayCecabundance.tif"))
clayCec <- calc(clayCec, fun = function(x){round(x * fixedeffs['scaleCLYPPT:scaleCECSOL'], digits = 2)},
              filename = file.path(savefolder, "clayCecabundance.tif"), overwrite = TRUE)

siltCec <- overlay(silt, cec, fun = createInteractionCoef,
                   filename = file.path(savefolder, "siltCecabundance.tif"))
siltCec <- calc(siltCec, fun = function(x){round(x * fixedeffs['scaleSLTPPT:scaleCECSOL'], digits = 2)},
                filename = file.path(savefolder, "siltCecabundance.tif"), overwrite = TRUE)

cecOrgC <- overlay(cec, orgC, fun = createInteractionCoef,
                 filename = file.path(savefolder, "cecCabundance.tif"))
ceOrgC <- calc(cecOrgC, fun = function(x){round(x * fixedeffs['scaleCECSOL:scaleORCDRC'], digits = 2)},
              filename = file.path(savefolder, "cecCabundancecoef.tif"), overwrite = TRUE)

print("Climate....")
bio1bio15 <- overlay(bio1, bio15, fun = createInteractionCoef, 
                  filename = file.path(savefolder, "bio1bio15abundance.tif"))
bio1bio15 <- calc(bio1bio15, fun = function(x){round(x * fixedeffs['bio10_1_scaled:bio10_15_scaled'], digits = 2)}, 
               filename = file.path(savefolder, "bio1bio15abundancecoef.tif"), overwrite = TRUE) 

bio1arid <- overlay(bio1, aridity, fun = createInteractionCoef, 
                     filename = file.path(savefolder, "bio1bio15abundance.tif"))
bio1arid <- calc(bio1arid, fun = function(x){round(x * fixedeffs['bio10_1_scaled:scaleAridity'], digits = 2)}, 
                  filename = file.path(savefolder, "bio1aridabundancecoef.tif"), overwrite = TRUE) 

bio1petsd <- overlay(bio1, petsd, fun = createInteractionCoef, 
                    filename = file.path(savefolder, "bio1petsdabundance.tif"))
bio1petsd <- calc(bio1petsd, fun = function(x){round(x * fixedeffs['bio10_1_scaled:ScalePETSD'], digits = 2)}, 
                 filename = file.path(savefolder, "bio1petsdabundancecoef.tif"), overwrite = TRUE) 

  
  
  
  
+ bio10_1_scaled:SnowMonths_cat +  
    
  bio10_15_scaled:SnowMonths_cat + SnowMonths_cat:ScalePETSD +  
  


print("Water retention....")
claybio15 <- overlay(bio15, clay, fun = createInteractionCoef, 
                     filename = file.path(savefolder, "claybio15abundance.tif"))
claybio15 <- calc(claybio15, fun = function(x){round(x * fixedeffs['scaleCLYPPT:bio10_15_scaled'], digits = 2)}, 
                  filename = file.path(savefolder, "claybio15abundancecoef.tif"), overwrite = TRUE) 

claypetsd <- overlay(clay, petsd, fun = createInteractionCoef, 
                     filename = file.path(savefolder, "petsdclayabundance.tif"))
claypetsd <- calc(claypetsd, fun = function(x){round(x * fixedeffs['scaleCLYPPT:ScalePETSD'], digits = 2)}, 
                  filename = file.path(savefolder, "petsdclayabundancecoef.tif"), overwrite = TRUE) 


#### Got here
## But didn't do the categorical variables

######## Add them all together
print("Adding all layers together")
cec <- calc(cec, fun = function(x) {round(x +  fixedeffs['(Intercept)'], digits = 2)}, 
            filename = file.path(savefolder, "CECrichnesscoef_plusInt.tif"))

f_together <- function(a, b, c, d, e, f, g, h, i, j){
  round(a + b +c +d +e +f +g + h+ i + j , digits = 2)
}

spR_finalraster <- overlay(bio1, bio4, bio12, bio15, ph, cec, orgC, 
                           phOrgC, bio1bio4, bio12bio15, fun = f_together, 
                           filename = file.path(savefolder, "spRFinalRasterMinusESA.tif"))
min(spR_finalraster)
max(spR_finalraster)
spR_finalraster
print("Done!") 