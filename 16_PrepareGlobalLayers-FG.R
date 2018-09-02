library(raster)
## To stop raster using the tmp dir which is slow and fills up
rasterOptions(tmpdir = "/work/phillips", chunksize = 524288, maxmemory = 134217728)

scaleGL <- function(layername, data, load = ".", save = "."){
  if(!(layername %in% names(data))) stop("Layername must be in data")
  
  print(layername)
  
  dat <- data[,which(names(data) == layername)]
  
  min <- min(dat, na.rm = TRUE)
  max <- max(dat, na.rm = TRUE)
  
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
  tif <- scale(tif)
  
  print("Saving tif")
  writeRaster(tif, save, format = "GTiff")
  # do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))
  
  
  rm(tif)
}



args <- commandArgs(trailingOnly = TRUE)

folder <- args[1] # data_dir
savefolder <- args[2] # output_dir
date <- args[3]
processed_dir <- args[4]

print(folder)
print(savefolder)
print(date)
print(processed_dir)

list.files(folder)


print("Loading dataset")

fg_richness <- read.csv(file = file.path(folder, paste("sites+FGRichness_", date, ".csv", sep = "")))



print("Calculating FG richness layers")

scaleGL(layername = "PHIHOX", data = fg_richness, load =  file.path(processed_dir, "PHIHOX_weighted.tif"),
        save = file.path(savefolder, 'PHIHOX_FGRichnessCutScaled.tif'))
scaleGL(layername = "ORCDRC", data = fg_richness, load =  file.path(processed_dir, "ORCDRC_weighted.tif"),
        save = file.path(savefolder, 'ORCDRC_FGRichnessCutScaled.tif'))
scaleGL(layername = "SLTPPT", data = fg_richness, load =  file.path(processed_dir, "SLTPPT_weighted.tif"),
        save = file.path(savefolder, 'SLTPPT_FGRichnessCutScaled.tif'))
scaleGL(layername = "CLYPPT", data = fg_richness, load =  file.path(processed_dir, "CLYPPT_weighted.tif"),
        save = file.path(savefolder, 'CLYPPT_FGRichnessCutScaled.tif'))
scaleGL(layername = "CECSOL", data = fg_richness, load =  file.path(processed_dir, "CECSOL_weighted.tif"),
      save = file.path(savefolder, 'CECSOL_FGRichnessCutScaled.tif'))

scaleGL(layername = 'bio10_1', data = fg_richness, load = file.path(processed_dir, 'CHELSA_bio10_1.tif'),
        save = file.path(savefolder, 'CHELSA_bio10_1_FGRichnessCutScaled.tif'))
scaleGL(layername = 'bio10_15', data = fg_richness, load = file.path(processed_dir, 'CHELSA_bio10_15.tif'),
        save = file.path(savefolder, 'CHELSA_bio10_15_FGRichnessCutScaled.tif'))
 
scaleGL(layername = 'Aridity', data = fg_richness, load = file.path(processed_dir, 'ai_yr_TIF.tif'),
         save = file.path(savefolder, 'Aridity_FGRichnessScaled.tif'))
scaleGL(layername = 'PET_SD', data = fg_richness, load = file.path(processed_dir, 'pet_he_SD.tif'),
         save = file.path(savefolder, 'PETSD_FGRichnessScaled.tif'))


print("Done!")

