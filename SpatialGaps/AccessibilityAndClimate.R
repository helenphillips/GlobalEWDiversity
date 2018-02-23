library(raster)
library(rgdal)

args <- commandArgs(trailingOnly = TRUE)

data_dir <- args[1]
output_dir <- args[2]

print(data_dir)
print(output_dir)

## To stop raster using the tmp dir which is slow and fills up
rasterOptions(tmpdir = "/work/phillips", chunksize = 524288, maxmemory = 134217728)

# GLdir <- "/data/idiv_sdiv"
# CLdir <- "/data/idiv_sdiv"
# outDir <- "/work/phillips"
## Load Accessibility layer
acc <- raster(file.path(data_dir, "accessibility_to_cities_2015_v1.tif"))
print("Loading accessibility layer")
## Load past climate layers
# precipitation
# p_precip <- raster(file.path(CLdir, "CHELSA_pr_1985-2010.tif"))
# Temperature
# p_temp <- raster(file.path(CLdir, "CHELSA_tas_1985-2010.tif"))
## Load future climate layers
# precipitation
f_precip <- raster(file.path(data_dir, "CHELSA_pr_1985-2050.tif"))
print("loaded precipitation layer")
# Temperature
f_temp <- raster(file.path(data_dir, "CHELSA_tas_1985-2050.tif"))
print("loaded temperature layer")


## They are different extents
f_temp <- crop(f_temp, acc, filename= file.path(output_dir, "CHELSA_tas_1985-2050_crop.tif"))
f_precip <- crop(f_precip, acc, filename= file.path(output_dir, "CHELSA_pr_1985-2050_crop.tif"))
print("cropped layers")

## If temperature change is above 20 (becuse I haven't divided by 10)
print("doing First subset")
acc[f_temp < 20] <- NA

print("doing second subset")
## and if rainfall is 
acc[f_precip > -20 & f_precip < 20] <- NA

## and if accessibility is very low
#  accessibility is in minutes. This is 24 hours
print("doing third subset")
acc[acc > 1440] <- NA


## make it into a mask
print("value to 1")
acc[!(is.na(acc))] <- 1

print("saving...")

### save it
rf <- writeRaster(acc,  filename=file.path(output_dir, "AccClimateMask.tif"), format="GTiff", overwrite=TRUE)