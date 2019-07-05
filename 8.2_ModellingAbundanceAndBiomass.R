args <- commandArgs(trailingOnly = TRUE)

data_out <- args[1] # output directory
data_in <- args[2] # data folder

# data_in <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\Results\\Revised"



##############

##### Functions
library(raster)
# require(biglm)
library(data.table)

rasterOptions(tmpdir = "/work/phillips", chunksize = 524288, maxmemory = 134217728)


# rasterOptions(tmpdir = "I:\\sDiv-PostDocs-Work\\Phillips\\tmp", chunksize = 524288, maxmemory = 134217728)



# regions <- c("africa", "asia", "europe", "latin_america", "north_america", "west_asia")

results <- file.path(data_in, "Biomass")
resultRaster <- "BiomassFinalRaster.tif"

print("Loading biomass rasters...")

africa <- as.vector(raster(file.path(results, "africa", resultRaster)))
asia <- as.vector(raster(file.path(results, "asia", resultRaster)))
europe <- as.vector(raster(file.path(results, "europe", resultRaster)))
latin_america <- as.vector(raster(file.path(results, "latin_america", resultRaster)))
north_america <- as.vector(raster(file.path(results, "north_america", resultRaster)))
west_asia <- as.vector(raster(file.path(results, "west_asia", resultRaster)))


biomass <- c(africa, asia, europe, latin_america, north_america, west_asia)

rm(list=c("africa", "asia", "europe", "latin_america", "north_america", "west_asia"))

############################################################
## ABUNDANCE
############################################################

results <- file.path(data_in, "Abundance")
resultRaster <- "AbundanceFinalRaster.tif"
print("Loading abundance rasters...")

africa <- as.vector(raster(file.path(results, "africa", resultRaster)))
asia <-  as.vector(raster(file.path(results, "asia", resultRaster)))
europe <- as.vector(raster(file.path(results, "europe", resultRaster)))
latin_america <- as.vector(raster(file.path(results, "latin_america", resultRaster)))
north_america <- as.vector(raster(file.path(results, "north_america", resultRaster)))
west_asia <- as.vector(raster(file.path(results, "west_asia", resultRaster)))

abundance <- c(africa, asia, europe, latin_america, north_america, west_asia)

rm(list=c("africa", "asia", "europe", "latin_america", "north_america", "west_asia"))

############################################################
## Plotting
############################################################
print("Creating data table")
div <- as.data.table(data.frame(abundance=abundance, biomass=biomass))


print("Removing the NAs")
div <- div[complete.cases(div),]

## Take a sample of 10%
print("Taking a sample of 10%")
ids <- sample(nrow(div), size = nrow(div)/10)
div <- div[ids,]

print("Saving data table")
save(div, file=file.path(data_out, "AbundanceBiomasDT.RData"))
print("Saving as CSV")
#write.csv(div,file=file.path(data_out, "AbundanceBiomasDT.csv"), row.names = FALSE)

fwrite(div, file.path(data_out, "AbundanceBiomasDT.csv"))

#load(file=file.path(data_out, "AbundanceBiomasDT.RData"))
print("Done!")
