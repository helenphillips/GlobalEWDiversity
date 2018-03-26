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
print("Loading accessibility layer")
acc <- raster(file.path(data_dir, "accessibility_to_cities_2015_v1.tif"))

## if accessibility is very low
#  accessibility is in minutes. This is 24 hours
print("doing subset based on acccesibility")
acc[acc > 720] <- NA



## Load future climate layers
# precipitation
print("loading precipitation layer")
f_precip <- raster(file.path(data_dir, "CHELSA_pr_1985-2050.tif"))

## They are different extents
print("cropping precipitation layer")
f_precip <- crop(f_precip, acc, filename= file.path(output_dir, "CHELSA_pr_1985-2050_crop.tif"), overwrite=TRUE)


# Temperature
print("loading temperature layer")
f_temp <- raster(file.path(data_dir, "CHELSA_tas_1985-2050.tif"))


## They are different extents
print("cropping temperature layer")
f_temp <- crop(f_temp, acc, filename= file.path(output_dir, "CHELSA_tas_1985-2050_crop.tif"), overwrite=TRUE)

## If temperature change is above 25 (becuse I haven't divided by 10)
## AND precitation is greater than abs. 50
print("doing subset based on precipitation and temp layer")
acc[(f_precip > -100 & f_precip < 100) & f_temp < 30] <- NA

# Removing precipitation layer
print("removing precipitation layer")
rm(f_precip)

# Removing temperature raster
print("removing temperature layer")
rm(f_temp)


## make it into a mask
print("value to 1")
acc[!(is.na(acc))] <- 1

print("saving...")

### save it
rf <- writeRaster(acc,  filename=file.path(output_dir, "AccClimateMask.tif"), format="GTiff", overwrite=TRUE)

print("Done!")