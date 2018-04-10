
########################################################
# 1. Load Libraries and Data
########################################################


library(raster)

########################################################
# 2. Set Working Directory or Cluster Info
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")

  climate_GLs <- "C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim"
  soil_GLs <- "I:\\sDiv-postdocs-work\\Phillips\\sWorm\\SpatialAnalysis\\SoilGrids\\1km"
  models <- "Models"
  }

if(Sys.info()["nodename"] == "DESKTOP"){
  setwd("D:\\Helens\\sWorm")
  
  climate_GLs <- "."
  soil_GLs <- "."
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
load(file.path(models, "richnessmodel_full.rds"))
load(file.path(models, "biomassmodel_full.rds"))
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
# 5. RICHNESS
#################################################
print("Creating richness raster")
print("Loading all rasters")
bio1 <- raster(file.path(GLs_folder, "CHELSA_bio10_1_RichnessCutScaled.tif"))
bio4 <- raster(file.path(GLs_folder, "CHELSA_bio10_4_RichnessCutScaled.tif"))
bio12 <- raster(file.path(GLs_folder, "CHELSA_bio10_12_RichnessCutScaled.tif"))
bio15 <- raster(file.path(GLs_folder, "CHELSA_bio10_15_RichnessCutScaled.tif"))

ph <- raster(file.path(GLs_folder,"PHIHOX_RichnessCutScaled.tif"))
# clay <- raster(file.path(soil_GLs,"CLYPPT_RichnessCutScaled.tif"))
cec <- raster(file.path(GLs_folder,"CECSOL_RichnessCutScaled.tif"))
orgC <- raster(file.path(GLs_folder,"ORCDRC_RichnessCutScaled.tif"))

esa <- raster(file.path(GLs_folder,"ESA_newValuesCropped.tif"))


fixedeffs <- summary(richness_model)$coefficients[,1]
# do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))


print("Calculating individual coefficients")
print("Climate coefficients....")
# fun <- function(x) { round(x * fixedeffs['bio10_1_scaled'], digits = 2 )}
# bio1 <- calc(bio1, fun, filename = file.path(savefolder, "bio10_1_richnesscoef.tif")) 
bio1 <- raster(file.path(GLs_folder, "bio10_1_richnesscoef.tif"))

#fun <- function(x) { round(x * fixedeffs['bio10_4_scaled'], digits = 2 )}
#bio4 <- calc(bio4, fun, filename = file.path(savefolder, "bio10_4_richnesscoef.tif")) 
bio4 <- raster(file.path(GLs_folder, "bio10_4_richnesscoef.tif"))


#fun <- function(x) { round(x * fixedeffs['bio10_12_scaled'], digits = 2 )}
#bio12 <- calc(bio12, fun, filename = file.path(savefolder, "bio10_12_richnesscoef.tif")) 
bio12 <- raster(file.path(GLs_folder, "bio10_12_richnesscoef.tif"))


#fun <- function(x) { round(x * fixedeffs['bio10_15_scaled'], digits = 2 )}
#bio15 <- calc(bio15, fun, filename = file.path(savefolder, "bio10_15_richnesscoef.tif")) 
bio15 <- raster(file.path(GLs_folder, "bio10_15_richnesscoef.tif"))

print("Soil coefficients....")
#fun <- function(x) { round(x * fixedeffs['scalePH'], digits = 2 ) }
#ph <- calc(ph, fun, filename = file.path(savefolder, "ph_richnesscoef.tif"))
ph<- raster(file.path(GLs_folder, "ph_richnesscoef.tif"))

#fun <- function(x) { round(x * fixedeffs['scaleCLYPPT'], digits = 2 ) }
#clay2 <- calc(clay, fun, filename = file.path(GLs_folder, "clay_richnesscoef.tif"))
# clay2 <- raster(file.path(soil_GLs, "clay_richnesscoef.tif"))

#fun <- function(x) { round(x * fixedeffs['scaleCECSOL'], digits = 2 ) }
#cec <- calc(cec, fun, filename = file.path(savefolder, "cec_richnesscoef.tif"))
cec <- raster(file.path(GLs_folder, "cec_richnesscoef.tif"))

#fun <- function(x) { round(x *fixedeffs['scaleORCDRC'], digits = 2 ) }
#orgC <- calc(orgC, fun, filename = file.path(savefolder, "orgc_richnesscoef.tif"))
orgC <- raster(file.path(GLs_folder, "orgc_richnesscoef.tif"))

print("Habitat cover...")
# esa[esa == 0] <- NA ## They have an NA pixel number 
# esa[esa == 60] <- 0 ## intercept is added later
# print("2 done")
# esa[esa == 50] <- round(fixedeffs['ESABroadleaf evergreen forest'], digits = 2)
# esa[esa == 70] <- round(fixedeffs['ESANeedleleaf evergreen forest'], digits = 2)
# print("4 done")
# esa[esa == 90] <- round(fixedeffs['ESAMixed forest'], digits = 2)
# esa[esa == 100] <- round(fixedeffs['ESATree open'], digits = 2)
# print("6 done")
# esa[esa == 110] <- round(fixedeffs['ESAHerbaceous with spare tree/shrub'], digits = 2)
# esa[esa == 120] <- round(fixedeffs['ESAShrub'], digits = 2)
# print("8 done")
# esa[esa == 130] <- round(fixedeffs['ESAHerbaceous'], digits = 2)
# esa[esa == 10] <- round(fixedeffs['ESAProduction - Herbaceous'], digits = 2)
# print("10 done. Last 3...")
# esa[esa == 12] <- round(fixedeffs['ESAProduction - Plantation'], digits = 2)
# esa[esa == 210] <- round(fixedeffs['ESAWater bodies'], digits = 2)
# do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))

## Put all habitat covers where we don't have the habitat cover in the model to NA
#esa[esa > 10] <- NA

print("Saving ESA layer")
#esa <- writeRaster(esa,  filename=file.path(savefolder, "ESA_coefs.tif"), format="GTiff", overwrite=TRUE)
#esa <- raster(file.path(GLs_folder, "ESA_coefs.tif"))
#esa
#str(esa)
#summary(esa)
#unique(getValues(esa))
print("Calculating interacting coefficients") 
print("Soil....Just 1")
# phOrgC <- overlay(ph, orgC, fun = createInteractionCoef, 
#                  filename = file.path(savefolder, "phCrichness.tif"))
# phOrgC <- calc(phOrgC, fun = function(x){round(x * fixedeffs['scalePH:scaleORCDRC'], digits = 2)}, 
#               filename = file.path(savefolder, "phCrichnesscoef.tif"), overwrite = TRUE) 
phOrgC <- raster(file.path(GLs_folder, "phCrichnesscoef.tif"))
  
print("Climate....")
#bio1bio4 <- overlay(bio1, bio4, fun = createInteractionCoef, 
#                  filename = file.path(savefolder, "bio1bio4richness.tif"))
#bio1bio4 <- calc(bio1bio4, fun = function(x){round(x * fixedeffs['bio10_1_scaled:bio10_4_scaled'], digits = 2)}, 
#               filename = file.path(savefolder, "bio1bio4richnesscoef.tif"), overwrite = TRUE) 
bio1bio4 <- raster(file.path(GLs_folder, "bio1bio4richnesscoef.tif"))


print("1 down, 1 to go")
#bio12bio15 <- overlay(bio12, bio15, fun = createInteractionCoef, 
#                    filename = file.path(savefolder, "bio12bio15richness.tif"))
#bio12bio15 <- calc(bio12bio15, fun = function(x){round(x * fixedeffs['bio10_12_scaled:bio10_15_scaled'], digits = 2)}, 
#                 filename = file.path(savefolder, "bio12bio15richnesscoef.tif"), overwrite = TRUE) 
bio12bio15 <- raster(file.path(GLs_folder, "bio12bio15richnesscoef.tif"))

# do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE))) 


######## ADd them all together
print("Cropping CHELSA layers")
bio1 <- crop(bio1, cec, filename= file.path(savefolder, "bio1richnesscoef_crop.tif"))
#bio1 <- raster(file.path(GLs_folder, "bio1richnesscoef_crop.tif"))

bio4 <- crop(bio4, cec, filename= file.path(savefolder, "bio4richnesscoef_crop.tif"))
#bio4 <- raster(file.path(GLs_folder, "bio4richnesscoef_crop.tif"))

bio12 <- crop(bio12, cec, filename= file.path(savefolder, "bio12richnesscoef_crop.tif"))
#bio12 <- raster(file.path(GLs_folder, "bio12richnesscoef_crop.tif"))

print("Three down, 3 to go")
bio15 <- crop(bio15, cec, filename= file.path(savefolder, "bio15richnesscoef_crop.tif"))
#bio15 <- raster(file.path(GLs_folder, "bio15richnesscoef_crop.tif"))

bio1bio4 <- crop(bio1bio4, cec,filename= file.path(savefolder, "bio1bio4richnesscoef_crop.tif"))
#bio1bio4 <- raster(file.path(GLs_folder, "bio1bio4richnesscoef_crop.tif"))

bio12bio15 <- crop(bio12bio15, cec, filename= file.path(savefolder, "bio12bio15richnesscoef_crop.tif"))
#bio12bio15 <- raster(file.path(GLs_folder, "bio12bio15richnesscoef_crop.tif"))
dim(bio12bio15)

print("Croping soil laters")
#phOrgC <- crop(phOrgC, esa, snap= "near", filename= file.path(savefolder, "phOrgCrichnesscoef_crop.tif"))
#ph <- crop(ph, esa,snap= "near", filename= file.path(savefolder, "phrichnesscoef_crop.tif"))
#cec <- crop(cec, esa, snap= "near",filename= file.path(savefolder, "cecrichnesscoef_crop.tif"))
#orgC <- crop(orgC, esa, snap= "near",filename= file.path(savefolder, "orgCrichnesscoef_crop.tif"))
ph
str(ph)
summary(ph)
dim(cec)
dim(ph)
# esa <- crop(esa, cec, filename= file.path(savefolder, "esarichnesscoef_crop.tif"))
#dim(esa)

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