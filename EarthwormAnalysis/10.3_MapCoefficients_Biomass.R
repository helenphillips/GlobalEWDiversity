
########################################################
# 1. Load Libraries and Data
########################################################


library(raster)
library(lme4)

########################################################
# 2. Set Working Directory or Cluster Info
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
  
  GLs_folder <- "I:\\sWorm\\ProcessedGLs\\Same_resolution\\regions"
  models <- "Models"
}else{ ## i.e. cluster
  args <- commandArgs(trailingOnly = TRUE)
  
  GLs_folder <- args[1] # GLs_dir
  models <- args[2] # models_dir
  savefolder <- args[3] # output_dir
  reg <- args[4] ## Which continent
  
  print(GLs_folder)
  print(models) 
  print(savefolder)
  print(reg)
  
  rasterOptions(tmpdir = "/work/phillips", chunksize = 524288, maxmemory = 134217728)
  
}

########################################################
# 2.5 Regions for analysis
########################################################

#################################################
# 3. Load in models
#################################################
print("Loading in the biodiversity models")
load(file.path(models, "biomassmodel_full_revised.rds"))


#################################################
# 4. Rerun model with different factor levels for ESA
#################################################
print("Re-running model with new ESA values....")
data <- biomass_model@frame
levels(data$ESA)[levels(data$ESA) == 'Broadleaf deciduous forest'] <- "60"
levels(data$ESA)[levels(data$ESA) == 'Broadleaf evergreen forest'] <- "50"
levels(data$ESA)[levels(data$ESA) == 'Needleleaf evergreen forest'] <- "70"
levels(data$ESA)[levels(data$ESA) == 'Mixed forest'] <- "90"
levels(data$ESA)[levels(data$ESA) == 'Herbaceous with spare tree/shrub'] <- "110"
levels(data$ESA)[levels(data$ESA) == 'Shrub'] <- "120"
levels(data$ESA)[levels(data$ESA) == 'Herbaceous'] <- "130"
levels(data$ESA)[levels(data$ESA) == 'Production - Herbaceous'] <- "10"
levels(data$ESA)[levels(data$ESA) == 'Production - Plantation'] <- "12"
# levels(data$ESA)[levels(data$ESA) == 'Cropland/Other vegetation mosaic'] <- "30"
levels(data$ESA)[levels(data$ESA) == "Bare area (unconsolidated"] <- "202"


mod <-  lmer(formula = biomass_model@call$formula, data = data, 
             control = lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))



####


if(!dir.exists(file.path(savefolder, reg))){
  dir.create(file.path(savefolder, reg))
}
#  data_out <- file.path(savefolder, reg)

#################################################
# 5. Biomass
#################################################
print("Creating biomass raster")
print("Loading all rasters")

print(file.path(GLs_folder, reg))

bio10_12_scaled <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_12_BiomassCutScaled.tif"))
dimensions <- dim(bio10_12_scaled)
resol <-res(bio10_12_scaled)
coordred <- crs(bio10_12_scaled)
exten <- extent(bio10_12_scaled)
bio10_12_scaled <- as.vector(bio10_12_scaled)

bio10_7_scaled <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_7_BiomassCutScaled.tif"))
bio10_7_scaled <- as.vector(bio10_7_scaled)


bio10_15_scaled <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_15_BiomassCutScaled.tif"))
bio10_15_scaled <- as.vector(bio10_15_scaled)


ScalePET <- raster(file.path(GLs_folder, reg, "PETyr_BiomassScaled.tif"))
ScalePET <- as.vector(ScalePET)


SnowMonths_cat <- raster(file.path(GLs_folder, reg, "Snow_newValues_WGS84.tif"))
SnowMonths_cat <- as.vector(SnowMonths_cat)
SnowMonths_cat <- as.factor(SnowMonths_cat)
levels(SnowMonths_cat)[levels(SnowMonths_cat) == "4"] <- "4plus"



scalePH <- raster(file.path(GLs_folder, reg,"PHIHOX_BiomassCutScaled.tif"))
scalePH <- as.vector(scalePH)

scaleCLYPPT <- raster(file.path(GLs_folder, reg, "CLYPPT_BiomassCutScaled.tif"))
scaleCLYPPT <- as.vector(scaleCLYPPT)

scaleSLTPPT <- raster(file.path(GLs_folder, reg, "SLTPPT_BiomassCutScaled.tif"))
scaleSLTPPT <- as.vector(scaleSLTPPT)

scaleCECSOL <- raster(file.path(GLs_folder, reg, "CECSOL_BiomassCutScaled.tif"))
scaleCECSOL <- as.vector(scaleCECSOL)

scaleORCDRC <- raster(file.path(GLs_folder, reg, "ORCDRC_BiomassCutScaled.tif"))
scaleORCDRC <- as.vector(scaleORCDRC)

ESA <- raster(file.path(GLs_folder, reg, "ESA_newValuesCropped.tif"))
ESA <- as.vector(ESA)
keep <- c(60, 50, 70, 90, 110, 120, 130, 10, 12, 202)
ESA <- ifelse(ESA %in% keep, ESA, NA)
ESA <- as.factor(ESA)

newdat <- data.frame(ESA = ESA,
                     scaleORCDRC = scaleORCDRC,
                     scaleCECSOL = scaleCECSOL,
                     scaleSLTPPT = scaleSLTPPT,
                     scaleCLYPPT = scaleCLYPPT,
                     scalePH = scalePH,
                     SnowMonths_cat = SnowMonths_cat,
                     ScalePET = ScalePET,
                     bio10_15_scaled = bio10_15_scaled,
                     bio10_12_scaled = bio10_12_scaled,
                     bio10_7_scaled = bio10_7_scaled)

rm(list=c("bio10_12_scaled", "bio10_15_scaled", "bio10_7_scaled", "SnowMonths_cat", "ScalePET",
          "scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL", "scaleORCDRC", "ESA"))
#############################################################
print("Splitting dataframe...")
library(data.table)
n <- 3000

letterwrap <- function(n, depth = 1) {
  args <- lapply(1:depth, FUN = function(x) return(LETTERS))
  x <- do.call(expand.grid, args = list(args, stringsAsFactors = F))
  x <- x[, rev(names(x)), drop = F]
  x <- do.call(paste0, x)
  if (n <= length(x)) return(x[1:n])
  return(c(x, letterwrap(n - length(x), depth = depth + 1)))
}

t <- nrow(newdat) %/% n 
alp <- letterwrap(t, depth = ceiling(log(t, base = 26)))
last <- alp[length(alp)]




print("1")
t <- rep(alp, each = n)
rm(alp)
# more <- letterwrap(1, depth = nchar(last) + 1)
more <- rep("Z", length = nchar(last) + 1)
implode <- function(..., sep='') {
  paste(..., collapse=sep)
}
more <- implode(more)

print("2")
newdat$z <- c(t, rep(more, times = (nrow(newdat) - length(t))))
rm(more)
rm(t)
rm(n)

print("3")
newdat_t = as.data.table(newdat)
rm(newdat)

gc()

print("4")
#system.time(
x <- split(newdat_t, f = newdat_t$z)
#)

rm(newdat_t)


print("Predicting values...")
# x <- split(newdat, (0:nrow(newdat) %/% 10000))  # modulo division


for(l in 1:length(x)){
  
  print(paste(l, "in", length(x), "iterations.."))
  
  res <- predict(mod, x[[l]], re.form = NA)
  write.table(res, file= file.path(savefolder, reg, "predictedValues.csv"),
              append=TRUE, row.names = FALSE,
              col.names = FALSE,
              sep = ',')
  
  x[[l]]$pred <- res
  
  write.table(x[[l]], file= file.path(savefolder, reg, "AllValues.csv"),
              append=TRUE, row.names = FALSE,
              col.names = FALSE,
              sep = ',')
  
}

res <- NULL
x <- NULL

#  length(res) == nrow(newdat)

# need number of rows of the original raster
# The resolution
# the extent
# the coord.ref
# dimensions
# resol
print("Loading csv of predicted values and converting to vector....")
predValues <- read.csv(file.path(savefolder, reg, "predictedValues.csv"), header = FALSE)
predValues <- as.vector(predValues$V1)

print("Converting to raster...")
print(dimensions[1])
print(dimensions[2])

print(dimensions[1] * dimensions[2])

# dimensions <- c(3032, 3074)
r <- matrix(predValues, nrow = dimensions[1], ncol = dimensions[2], byrow = TRUE)
r <- raster(r)

print("Adding in the raster information")

extent(r) <- exten
# ... and assign a projection
projection(r) <- coordred


# Save raster
print("Saving raster...")
r <- writeRaster(r,  filename=file.path(savefolder, reg, "BiomassFinalRaster.tif"), format="GTiff", overwrite=TRUE)


print("Done!")
