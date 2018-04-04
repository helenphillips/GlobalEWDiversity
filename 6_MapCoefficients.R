########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")

  climate_GLs <- "C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim"
  soil_GLs <- "I:\\sDiv-postdocs-work\\Phillips\\sWorm\\SpatialAnalysis\\SoilGrids\\1km"
  
  }

if(Sys.info()["nodename"] == "DESKTOP"){
  setwd("D:\\Helens\\sWorm")
  
  climate_GLs <- "."
  soil_GLs <- "."
}

########################################################
# 1. Load Libraries and Data
########################################################
library(raster)

source("Functions/FormatData.R")




#################################################
# 2. Loading in variables
#################################################

data_in <-"4_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
# loadin <- files[grep(date, files)]

if(!dir.exists("6_Data")){
  dir.create("6_Data")
}
dataOut <- "6_Data"


#################################################
# 3. Load in data
#################################################

richness <- read.csv(file.path(data_in, paste("sitesRichness_",date,".csv", sep = "")))
biomass <- read.csv(file.path(data_in, paste("sitesBiomass_",date,".csv", sep = "")))
abundance <- read.csv(file.path(data_in, paste("sitesAbundance_",date,".csv", sep = "")))

richness <- droplevels(SiteLevels(richness))
biomass <- droplevels(SiteLevels(biomass))
abundance <- droplevels(SiteLevels(abundance))



#################################################
# 4. Load in models
#################################################

models <- "Models"


load(file.path(models, "richnessmodel_full.rds"))
load(file.path(models, "biomassmodel_full.rds"))
load(file.path(models, "abundancemodel_full.rds"))

#################################################
# 5. RICHNESS
#################################################

bio15 <- raster(file.path(climate_GLs, "CHELSA_bio10_15_RichnessCutScaled.tif"))

ph <- raster(file.path(soil_GLs,"PHIHOX_RichnessCutScaled.tif"))
clay <- raster(file.path(soil_GLs,"CLYPPT_RichnessCutScaled.tif"))
cec <- raster(file.path(soil_GLs,"CECSOL_RichnessCutScaled.tif"))
orgC <- raster(file.path(soil_GLs,"ORCDRC_RichnessCutScaled.tif"))


fixedeffs <- summary(richness_model)$coefficients[,1]
do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))

createCoef <- function(X, coef, ...){
  t <- x * coef
  round(t, digits = 2)
}

fun <- function(x) { round(x * fixedeffs['bio10_15_scaled'], digits = 2 )}
bio15 <- calc(bio15, fun, filename = file.path(soil_GLs, "bio10_15_richnesscoef.tif")) 
# bio15 <- raster(file.path(climate_GLs, "bio10_15_richnesscoef.tif"))

fun <- function(x) { round(x * fixedeffs['scalePH'], digits = 2 ) }
ph2 <- calc(ph, fun, filename = file.path(soil_GLs, "ph_richnesscoef.tif"))
# ph2 <- raster(file.path(soil_GLs, "ph_richnesscoef.tif"))

fun <- function(x) { round(x * fixedeffs['scaleCLYPPT'], digits = 2 ) }
clay2 <- calc(clay, fun, filename = file.path(soil_GLs, "clay_richnesscoef.tif"))
# clay2 <- raster(file.path(soil_GLs, "clay_richnesscoef.tif"))

fun <- function(x) { round(x * fixedeffs['scaleCECSOL'], digits = 2 ) }
cec2 <- calc(cec, fun, filename = file.path(soil_GLs, "cec_richnesscoef.tif"))
#cec2 <- calc(file.path(soil_GLs, "cec_richnesscoef.tif"))

fun <- function(x) { round(x *fixedeffs['scaleORCDRC'], digits = 2 ) }
orgC2 <- calc(orgC, fun, filename = file.path(soil_GLs, "orgc_richnesscoef.tif"))
#orgC2 <- raster(file.path(soil_GLs, "orgc_richnesscoef.tif"))


do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))

createInteractionCoef <- function(x, y){
  round(x * y, digits = 2)
}
 
phOrgC <- overlay(ph, orgC, fun = createInteractionCoef, 
                  filename = file.path(soil_GLs, "phCrichnesscoef.tif"))
phOrgC <- calc(phOrgC, fun = function(x){round(x * fixedeffs['scalePH:scaleORCDRC'], digits = 2)}, 
               filename = file.path(soil_GLs, "phCrichnesscoef_f.tif"), overwrite = TRUE) 
# phOrgC <- raster(file.path(soil_GLs, "phCrichnesscoef.tif"))
  

do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE))) 

claycec <- overlay(clay, cec, fun = createInteractionCoef, 
                     filename = file.path(soil_GLs, "cecclayrichnesscoef.tif"))
  # claycec <- raster(file.path(soil_GLs, "cecclayrichnesscoef.tif"))
do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))
claycec <- calc(claycec, fun = function(x){round(x * fixedeffs['scaleCLYPPT:scaleCECSOL'], digits = 2)}, 
               filename = file.path(soil_GLs, "cecclayrichnesscoef_f.tif"), overwrite = TRUE) 

######## ADd them all together
  
bio15 <- crop(bio15, clay2, filename= file.path(soil_GLs, "bio15richnesscoef_crop.tif"))
  
bio15 <- calc(bio15, fun = function(x) {round(x +  fixedeffs['(Intercept)'], digits = 2)}, 
              filename = file.path(soil_GLs, "bio15richnesscoef_cropInt.tif"))
 f_together <- function(a, b, c, d, e, f, g){
    round(a + b +c +d +e +f +g, digits = 2)
  }
  
spR_finalraster <- overlay(bio15, clay2, ph2, cec2, orgC2, phOrgC, claycec, fun = f_together, 
                             filename = file.path(soil_GLs, "spRFinalRaster.tif"))
 