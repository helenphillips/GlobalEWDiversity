
library(raster)

scaleGL <- function(layername, data, load = ".", save = "."){
  if(!(layername %in% names(data))) stop("Layername must be in data")
  
  dat <- data[,which(names(data) == layername)]
  
  min <- min(dat, na.rm = TRUE)
  max <- max(dat, na.rm = TRUE)
  
  tif <- raster(load)
  
  dividby10 <- c("bio10_1","bio10_5", 
                 "bio10_6", "bio10_7","bio10_8","bio10_9","bio10_10",
                 "bio10_11", "PHIHOX", "ORCDRC")
  
  
  if(layername %in% dividby10){
    tif <- tif/10
  }
  
  tif[tif<min] <- min
  tif[tif>max] <- max
  
  
  tif <- scale(tif)
  
  writeRaster(tif, save, format = "GTiff")
  do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))
  

  rm(tif)
}


calcualteWeightedMeanGL <- function(folder, savefolder){
  tifs <- c("PHIHOX", "CLYPPT", "SLTPPT", "CECSOL", "ORCDRC") ## "SNDPPT",  "TAXNWRB_1"
  layers <- c("sl1","sl2", "sl3", "sl4")
  weight <- c(0.001, 0.05, 0.1, 0.15)
  filenames <- list.files(folder)
  
  for(t in tifs){
    dl <- filenames[grep(t, filenames)]
    dl <- dl[grep(paste(layers,collapse="|"), dl)]
    
    
    s <- raster::stack(file.path(folder, dl))
    ## writeRaster(s, "hdf8_EVI.TIF")
    
    weighted.mean(x = s, w = weight, na.rm=FALSE,
                  filename=file.path(savefolder, paste(t, "_weighted.tif", sep ="")), format = "GTiff")
    do.call(file.remove, list(list.files(dirname(rasterTmpFile()), full.names = TRUE)))
  }
  
}

folder <- 'D:/Helens/sWorm'
savefolder <- folder

richness <- read.csv(file.path(folder, 'sitesRichness_2017-12-04.csv'))
abundance <- read.csv(file.path(folder, 'sitesAbundance_2017-12-04.csv'))
biomass <- read.csv()


calcualteWeightedMeanGL(folder, savefolder)

## Richness

scaleGL(layername = 'bio10_15', data = richness, load = file.path(folder, 'CHELSA_bio10_15.tif'), 
        save = file.path(folder, 'CHELSA_bio10_15_RichnessCutScaled.tif'))


scaleGL(layername = "ORCDRC", data = richness, load =  file.path(folder, "ORCDRC_weighted.tif"), 
        save = file.path(folder, 'ORCDRC_RichnessCutScaled.tif'))
scaleGL(layername = "PHIHOX", data = richness, load =  file.path(folder, "PHIHOX_weighted.tif"), 
        save = file.path(folder, 'PHIHOX_RichnessCutScaled.tif'))
scaleGL(layername = "CLYPPT", data = richness, load =  file.path(folder, "CLYPPT_weighted.tif"), 
        save = file.path(folder, 'CLYPPT_RichnessCutScaled.tif'))
scaleGL(layername = "SLTPPT", data = richness, load =  file.path(folder, "SLTPPT_weighted.tif"), 
        save = file.path(folder, 'SLTPPT_RichnessCutScaled.tif'))
scaleGL(layername = "CECSOL", data = richness, load =  file.path(folder, "CECSOL_weighted.tif"), 
        save = file.path(folder, 'CECSOL_RichnessCutScaled.tif'))
## Abundance

scaleGL(layername = 'bio10_10', data = abundance, load = file.path(folder, 'CHELSA_bio10_10.tif'), 
        save = file.path(folder, 'CHELSA_bio10_10_RichnessCutScaled.tif'))
scaleGL(layername = 'bio10_18', data = abundance, load = file.path(folder, 'CHELSA_bio10_18.tif'), 
        save = file.path(folder, 'CHELSA_bio10_18_RichnessCutScaled.tif'))





  bio10_8
bio10_2
bio10_15
bio10_5
bio10_10
bio10_18
  
  
PHIHOX
CLYPPT
SLTPPT
CECSOL
ORCDRC

