library(spatial.tools)
load <- file.path("C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim", 
          paste("CHELSA_bio10_18.tif", sep=""))


scaleGL <- function(layername , data, load = ".", save = "."){
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
  
  tif <- ifelse(tif < min, min, tif)
  tif <- ifelse(tif > max, max, tif)
  
  tif <- scale(tif)
  
  writeRaster(tif, save, format = "GTiff")
  rm(tif)
}


rst_test <- function(rast, min) 
{ 
  return(ifelse(rast < min, min, rast)) 
} 

# Fill in the rasters you intend to use here: 
raster_check <- 
  rasterEngine(rast=tif, fun=rst_test, args = list(min=min)) 


PHIHOX
CLYPPT
SLTPPT
CECSOL
ORCDRC

bio10_8
bio10_2
bio10_15
bio10_5
bio10_10
bio10_18