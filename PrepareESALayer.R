library(raster)


args <- commandArgs(trailingOnly = TRUE)

GLs_folder <- args[1] # GLs_dir
## Something like "/data/idiv_sdiv/sworm/ProccessGLs"
savefolder <- args[2] # output_dir
## Something like "/data/idiv_sdiv/sworm/ProccessGLs"

print(GLs_folder)
print(savefolder)

rasterOptions(tmpdir = "/work/phillips", chunksize = 524288, maxmemory = 134217728)



print("Loading ESA raster")
r <- raster(file.path(GLs_folder, "ESACCI-LC-L4-LCCS-Map-1000m-P5Y-2010-v1.6.1_NEW.tif"))
# getValues(r)
# plot(r)

print("Cropping based on soil layer")
ph <- raster(file.path(GLs_folder,"PHIHOX_RichnessCutScaled.tif"))
r <- crop(r, ph, filename= file.path(savefolder, "ESA_cropped.tif"), overwrite=TRUE)

print("Changing values")
r[r == 11] <- 10 # Production herbs
r[r == 40] <- 30
print("2 done")
r[r == 61] <- 60
r[r == 62] <- 60
print("4 done")
r[r == 71] <- 70
r[r == 72] <- 70
print("6 done")
r[r == 81] <- 80
r[r == 82] <- 80
print("8 done")
r[r == 121] <- 120
r[r == 122] <- 120
print("10 done")
r[r == 151] <- 150
r[r == 152] <- 150
r[r == 153] <- 150
print("12 done")
r[r == 202] <- 200




print("Saving!")
r <- writeRaster(r,  filename=file.path(savefolder, "ESA_newValuesCropped.tif"), format="GTiff", overwrite=TRUE)

print("Done!!!")