
library(raster)
## To stop raster using the tmp dir which is slow and fills up
rasterOptions(tmpdir = "/work/phillips", chunksize = 524288, maxmemory = 134217728)

scaleGL <- function(layername, data, load = ".", save = "."){
  if(!(layername %in% names(data))) stop("Layername must be in data")
  
  print(layername)
  
  dat <- data[,which(names(data) == layername)]
  
  min <- min(dat, na.rm = TRUE)
  max <- max(dat, na.rm = TRUE)
  
  mean_val <- mean(dat, na.rm = TRUE)
  sd_val <- sd(dat, na.rm = TRUE)

  tif <- raster(load)
  
  dividby10 <- c("bio10_1","bio10_5", 
                 "bio10_6", "bio10_7","bio10_8","bio10_9","bio10_10",
                 "bio10_11", "PHIHOX", "ORCDRC")
  
  
  if(layername %in% dividby10){
    print("Dividing by 10")
    tif <- tif/10
  }
  
  print("Capping minimum value")
  tif[tif<min] <- min
  print("Capping maximum value")
  tif[tif>max] <- max
  


  print("Scaling tif")
  #fun <- function(x, mean_val) { x - mean_val }
  #tif <- calc(tif, fun)
  tif <- tif - mean_val
  
  #fun <- function(x, sd_val) { x / sd_val }
  #tif <- calc(tif, fun)
  tif <- tif / sd_val

  
  
  
  print("Saving tif")
  writeRaster(tif, save, format = "GTiff")
  # do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))
  

  rm(tif)
}


calculateWeightedMeanGL <- function(folder, savefolder){
  # tifs <- c("PHIHOX", "CLYPPT", "SLTPPT", "CECSOL", "ORCDRC") ## "SNDPPT",  "TAXNWRB_1"
  tifs <- "ORCDRC"
  print(tifs)
  layers <- c("sl1","sl2", "sl3", "sl4")
  weight <- c(0.001, 0.05, 0.1, 0.15)
  filenames <- list.files(folder)
  
  for(t in tifs){
    print(paste("Calculating weighted mean of", t))
    dl <- filenames[grep(t, filenames)]
    dl <- dl[grep(paste(layers,collapse="|"), dl)]
    
    print(paste("Creating stack of", file.path(folder, dl)))
    s <- raster::stack(file.path(folder, dl))
    ## writeRaster(s, "hdf8_EVI.TIF")
    
    print("Doing calculations")
    weighted.mean(x = s, w = weight, na.rm=FALSE,
                  filename=file.path(savefolder, paste(t, "_weighted.tif", sep ="")), format = "GTiff")
    rm(s)
   #  do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))
  }
  
}


args <- commandArgs(trailingOnly = TRUE)

folder <- args[1] # data_dir
savefolder <- args[2] # output_dir
date <- args[3]
site_dir <- args[4]

print(folder)
print(savefolder)
print(date)
print(site_dir)

list.files(folder)

# folder <- 'D:/Helens/sWorm'
# savefolder <- folder
# date <- date

# print("Reclassifying Snow Layer")
# 
# r <- raster(file.path(folder, "snow_2015_sum.tif"))
# 
# print("Changing values")
# r[r > 3] <- 4 # Anything above 3
# 
# print("Saving!")
# r <- writeRaster(r,  filename=file.path(savefolder, "Snow_newValues.tif"), format="GTiff", overwrite=TRUE)
# 
# 
# 

print("Loading datasets")
richness <- read.csv(file.path(site_dir, paste('sitesRichness_', date, '.csv', sep="")))
# abundance <- read.csv(file.path(site_dir, paste('sitesAbundance_', date, '.csv', sep="")))
# biomass <- read.csv(file.path(site_dir, paste('sitesBiomass_', date, '.csv', sep="")))

# print("Calculating the weighted mean of soil data")
# calculateWeightedMeanGL(folder, savefolder)

 print("Calculating richness layers")
# Richness

# 
#   scaleGL(layername = "ORCDRC", data = richness, load =  file.path(folder, "ORCDRC_weighted.tif"),
#          save = file.path(savefolder, 'ORCDRC_RichnessCutScaled.tif'))
#   scaleGL(layername = "PHIHOX", data = richness, load =  file.path(folder, "PHIHOX_weighted.tif"),
#          save = file.path(savefolder, 'PHIHOX_RichnessCutScaled.tif'))
#   scaleGL(layername = "CLYPPT", data = richness, load =  file.path(folder, "CLYPPT_weighted.tif"),
#          save = file.path(savefolder, 'CLYPPT_RichnessCutScaled.tif'))
#   scaleGL(layername = "SLTPPT", data = richness, load =  file.path(folder, "SLTPPT_weighted.tif"),
#          save = file.path(savefolder, 'SLTPPT_RichnessCutScaled.tif'))
#   scaleGL(layername = "CECSOL", data = richness, load =  file.path(folder, "CECSOL_weighted.tif"),
#          save = file.path(savefolder, 'CECSOL_RichnessCutScaled.tif'))
#  
  # scaleGL(layername = 'bio10_1', data = richness, load = file.path(folder, 'CHELSA_bio10_1.tif'),
  #        save = file.path(savefolder, 'CHELSA_bio10_1_RichnessCutScaled.tif'))
  # scaleGL(layername = 'bio10_4', data = richness, load = file.path(folder, 'CHELSA_bio10_4.tif'),
  #        save = file.path(savefolder, 'CHELSA_bio10_4_RichnessCutScaled.tif'))
  # scaleGL(layername = 'bio10_7', data = richness, load = file.path(folder, 'CHELSA_bio10_7.tif'),
  #        save = file.path(savefolder, 'CHELSA_bio10_7_RichnessCutScaled.tif'))
  # 
   scaleGL(layername = 'bio10_12', data = richness, load = file.path(folder, 'CHELSA_bio10_12.tif'),
         save = file.path(savefolder, 'CHELSA_bio10_12_RichnessCutScaled.tif'))
  scaleGL(layername = 'bio10_15', data = richness, load = file.path(folder, 'CHELSA_bio10_15.tif'),
         save = file.path(savefolder, 'CHELSA_bio10_15_RichnessCutScaled.tif'))
 
  scaleGL(layername = 'Aridity', data = richness, load = file.path(folder, 'ai_yr_TIF.tif'),
          save = file.path(savefolder, 'Aridity_RichnessScaled.tif'))
  scaleGL(layername = 'PETyr', data = richness, load = file.path(folder, 'pet_he_yr_TIF.tif'),
          save = file.path(savefolder, 'PETyr_RichnessScaled.tif'))
  scaleGL(layername = 'PET_SD', data = richness, load = file.path(folder, 'pet_he_SD.tif'),
          save = file.path(savefolder, 'PETSD_RichnessScaled.tif'))
  scaleGL(layername = 'elevation', data = richness, load = file.path(folder, 'elevation.tif'),
          save = file.path(savefolder, 'elevation_RichnessScaled.tif'))
 
 ## Abundance
# print("Calculating abundance layers")
#   scaleGL(layername = "ORCDRC", data = abundance, load =  file.path(folder, "ORCDRC_weighted.tif"), 
#          save = file.path(savefolder, 'ORCDRC_AbundanceCutScaled.tif'))
#   scaleGL(layername = "PHIHOX", data = abundance, load =  file.path(folder, "PHIHOX_weighted.tif"), 
#          save = file.path(savefolder, 'PHIHOX_AbundanceCutScaled.tif'))
#   scaleGL(layername = "CLYPPT", data = abundance, load =  file.path(folder, "CLYPPT_weighted.tif"), 
#          save = file.path(savefolder, 'CLYPPT_AbundanceCutScaled.tif'))
#   scaleGL(layername = "SLTPPT", data = abundance, load =  file.path(folder, "SLTPPT_weighted.tif"), 
#          save = file.path(savefolder, 'SLTPPT_AbundanceCutScaled.tif'))
#   scaleGL(layername = "CECSOL", data = abundance, load =  file.path(folder, "CECSOL_weighted.tif"), 
#          save = file.path(savefolder, 'CECSOL_AbundanceCutScaled.tif'))
 
#    scaleGL(layername = 'bio10_1', data = abundance, load = file.path(folder, 'CHELSA_bio10_1.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_1_AbundanceCutScaled.tif'))
#   scaleGL(layername = 'bio10_4', data = abundance, load = file.path(folder, 'CHELSA_bio10_4.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_4_AbundanceCutScaled.tif'))
#   scaleGL(layername = 'bio10_7', data = abundance, load = file.path(folder, 'CHELSA_bio10_7.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_7_AbundanceCutScaled.tif'))
#   scaleGL(layername = 'bio10_12', data = abundance, load = file.path(folder, 'CHELSA_bio10_12.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_12_AbundanceCutScaled.tif'))
#  scaleGL(layername = 'bio10_15', data = abundance, load = file.path(folder, 'CHELSA_bio10_15.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_15_AbundanceCutScaled.tif'))
 
#  scaleGL(layername = 'Aridity', data = abundance, load = file.path(folder, 'ai_yr_TIF.tif'),
#          save = file.path(savefolder, 'Aridity_AbundanceScaled.tif'))
#  scaleGL(layername = 'PETyr', data = abundance, load = file.path(folder, 'pet_he_yr_TIF.tif'),
#          save = file.path(savefolder, 'PETyr_AbundanceScaled.tif'))
#  scaleGL(layername = 'PET_SD', data = abundance, load = file.path(folder, 'pet_he_SD.tif'),
#          save = file.path(savefolder, 'PETSD_AbundanceScaled.tif'))
# scaleGL(layername = 'elevation', data = abundance, load = file.path(folder, 'elevation.tif'),
#         save = file.path(savefolder, 'elevation_AbundanceScaled.tif'))

# # Biomass
# print("Calculating biomass layers")
#  scaleGL(layername = "ORCDRC", data = biomass, load =  file.path(folder, "ORCDRC_weighted.tif"), 
#          save = file.path(savefolder, 'ORCDRC_BiomassCutScaled.tif'))
#  scaleGL(layername = "PHIHOX", data = biomass, load =  file.path(folder, "PHIHOX_weighted.tif"), 
#          save = file.path(savefolder, 'PHIHOX_BiomassCutScaled.tif'))
#  scaleGL(layername = "CLYPPT", data = biomass, load =  file.path(folder, "CLYPPT_weighted.tif"), 
#          save = file.path(savefolder, 'CLYPPT_BiomassCutScaled.tif'))
#  scaleGL(layername = "SLTPPT", data = biomass, load =  file.path(folder, "SLTPPT_weighted.tif"), 
#          save = file.path(savefolder, 'SLTPPT_BiomassCutScaled.tif'))
#  scaleGL(layername = "CECSOL", data = biomass, load =  file.path(folder, "CECSOL_weighted.tif"), 
#          save = file.path(savefolder, 'CECSOL_BiomassCutScaled.tif'))
 
#  scaleGL(layername = 'bio10_1', data = biomass, load = file.path(folder, 'CHELSA_bio10_1.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_1_BiomassCutScaled.tif'))
#  scaleGL(layername = 'bio10_4', data = biomass, load = file.path(folder, 'CHELSA_bio10_4.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_4_BiomassCutScaled.tif'))
#  scaleGL(layername = 'bio10_7', data = biomass, load = file.path(folder, 'CHELSA_bio10_7.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_7_BiomassCutScaled.tif'))
#  scaleGL(layername = 'bio10_12', data = biomass, load = file.path(folder, 'CHELSA_bio10_12.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_12_BiomassCutScaled.tif'))
#  scaleGL(layername = 'bio10_15', data = biomass, load = file.path(folder, 'CHELSA_bio10_15.tif'), 
#          save = file.path(savefolder, 'CHELSA_bio10_15_BiomassCutScaled.tif'))
 
#  scaleGL(layername = 'Aridity', data = biomass, load = file.path(folder, 'ai_yr_TIF.tif'),
#          save = file.path(savefolder, 'Aridity_BiomassScaled.tif'))
#  scaleGL(layername = 'PETyr', data = biomass, load = file.path(folder, 'pet_he_yr_TIF.tif'),
#          save = file.path(savefolder, 'PETyr_BiomassScaled.tif'))
#  scaleGL(layername = 'PET_SD', data = biomass, load = file.path(folder, 'pet_he_SD.tif'),
#          save = file.path(savefolder, 'PETSD_BiomassScaled.tif'))
# scaleGL(layername = 'elevation', data = biomass, load = file.path(folder, 'elevation.tif'),
#         save = file.path(savefolder, 'elevation_BiomassScaled.tif'))

print("Done!")

