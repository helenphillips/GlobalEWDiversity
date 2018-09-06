
args <- commandArgs(trailingOnly = TRUE)

datafolder <- args[1] # data_in
savefolder <- args[2] # output_dir


print(datafolder)
print(savefolder)


library(raster)
rasterOptions(tmpdir = "/work/phillips", chunksize = 524288, maxmemory = 134217728)



##################################
# Create output document
##################################

print("Creating an output document")

filename <- file.path(savefolder, "MapsAnalysed.txt")
cat("Top and bottom 20 percent", file = filename)
cat("\n\n", file = filename, append = TRUE)

##################################
# Process the rasters
##################################

print("Processing the rasters....")


for(resultRaster in c("spRFinalRaster.tif", "BiomassFinalRaster.tif"))
{

  print(resultRaster)
  cat(resultRaster, file = filename, append = TRUE)
  cat("\n\n", file = filename, append = TRUE)
  
  if(resultRaster == "spRFinalRaster.tif"){results <- file.path(datafolder, "Richness")}
  if(resultRaster == "BiomassFinalRaster.tif"){results <- file.path(datafolder, "Biomass")}
  
  print(results)
  
  africa <- raster(file.path(results, "africa", resultRaster))
  if(resultRaster == "spRFinalRaster.tif"){
    africa <- exp(africa)
  } else { africa <- exp(africa) - 1}
  
  asia <-  raster(file.path(results, "asia", resultRaster))
  if(resultRaster == "spRFinalRaster.tif"){
    asia <- exp(asia)
  } else { asia <- exp(asia) - 1}
  
  europe <- raster(file.path(results, "europe", resultRaster))
  if(resultRaster == "spRFinalRaster.tif"){
    europe <- exp(europe)
  } else { europe <- exp(europe) - 1}
  
  latin_america <- raster(file.path(results, "latin_america", resultRaster))
  if(resultRaster == "spRFinalRaster.tif"){
    latin_america <- exp(latin_america)
  } else { latin_america <- exp(latin_america) - 1}
  
  north_america <- raster(file.path(results, "north_america", resultRaster))
  if(resultRaster == "spRFinalRaster.tif"){
    north_america <- exp(north_america)
  } else { north_america <- exp(north_america) - 1}
  
  west_asia <- raster(file.path(results, "west_asia", resultRaster))
  if(resultRaster == "spRFinalRaster.tif"){
    west_asia <- exp(west_asia)
  } else { west_asia <- exp(west_asia) - 1}
  

  minV <-min(c(minValue(africa),
             minValue(asia),
             minValue(europe),
             minValue(latin_america),
             minValue(north_america),
             minValue(west_asia)))
  
  maxV <-max(c(maxValue(africa),
             maxValue(asia),
             maxValue(europe),
             maxValue(latin_america),
             maxValue(north_america),
             maxValue(west_asia)))


  #####
  # Calculating the top X percent
  #####
  
  top20 <- maxV - (maxV * 0.2)
  bottom20 <- minV + (maxV * 0.2)
  
  regions_all <- c("africa", "asia", "europe", "latin_america", "north_america", "west_asia")
  
  dat <- as.data.frame(matrix(NA, nrow = length(regions_all), ncol = 4))
  dat[,1] <- regions_all
  names(dat) <- c("Area", "TotalCells_N", "TopCells_N", "BottomsCells_N")
  for(reg in 1:length(regions_all)){
    
    print(regions_all[reg])
    
    print("reclassifying bottom")
    m <- c(minV, bottom20, 1,  bottom20, maxV, 0)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    region <- reclassify(get(regions_all[reg]), rclmat)  ## "get" converts a string to a object name
    
    print("summing bottom")
    dat$BottomsCells_N[reg] <- sum(region[], na.rm = TRUE)
    
    print("reclassifying top")
    m <- c(minV, top20, 0,  top20, maxV, 1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    region <- reclassify(get(regions_all[reg]), rclmat)  ## "get" converts a string to a object name
    
    print("summing top")          
    dat$TopCells_N[reg] <- sum(region[], na.rm = TRUE)
    
    print("reclassifying total area") 
    m <- c(minV, maxV, 1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    region <- reclassify(get(regions_all[reg]), rclmat)
    print("summing total area")
    dat$TotalCells_N[reg] <- sum(region[], na.rm = TRUE)
   
  }
  
  cat(paste(resultRaster, "\n"), file = filename, append = TRUE)
  capture.output(dat, file = filename, append = TRUE)
  
}

print("Done!")
