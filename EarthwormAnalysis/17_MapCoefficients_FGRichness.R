
########################################################
# 1. Load Libraries and Data
########################################################


library(raster)

########################################################
# 2. Set Working Directory or Cluster Info
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
  
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


#################################################
# 3. Load in models
#################################################
print("Loading in the biodiversity models")
load(file.path(models, "fgrichnessmodel.rds"))



if(!dir.exists(file.path(savefolder, reg))){
  dir.create(file.path(savefolder, reg))
}
#  data_out <- file.path(savefolder, reg)



#################################################
# 4. Rerun model with different factor levels for ESA
#################################################
print("Re-running model with new ESA values....")
data <- fgrichness_model@frame
levels(data$ESA)[levels(data$ESA) == 'Broadleaf deciduous forest'] <- "60"
levels(data$ESA)[levels(data$ESA) == 'Broadleaf evergreen forest'] <- "50"
levels(data$ESA)[levels(data$ESA) == 'Needleleaf evergreen forest'] <- "70"
levels(data$ESA)[levels(data$ESA) == 'Mixed forest'] <- "90"
levels(data$ESA)[levels(data$ESA) == 'Herbaceous with spare tree/shrub'] <- "110"
levels(data$ESA)[levels(data$ESA) == 'Shrub'] <- "120"
levels(data$ESA)[levels(data$ESA) == 'Herbaceous'] <- "130"
levels(data$ESA)[levels(data$ESA) == 'Production - Herbaceous'] <- "10"
levels(data$ESA)[levels(data$ESA) == 'Production - Plantation'] <- "12"
levels(data$ESA)[levels(data$ESA) == 'Cropland/Other vegetation mosaic'] <- "30"


mod <-  glmer(formula = fgrichness_model@call$formula, data = data, family = "poisson",
              control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))


#################################################
# 5. FG Richness
#################################################
print("Creating Functional group richness raster")
print("Loading all rasters")

print(file.path(GLs_folder, reg))




bio10_1_scaled <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_1_FGRichnessCutScaled.tif"))
dimensions <- dim(bio10_1_scaled)
resol <-res(bio10_1_scaled)
coordred <- crs(bio10_1_scaled)
exten <- extent(bio10_1_scaled)
bio10_1_scaled <- as.vector(bio10_1_scaled)


bio10_15_scaled <- raster(file.path(GLs_folder, reg, "CHELSA_bio10_15_FGRichnessCutScaled.tif"))
bio10_15_scaled <- as.vector(bio10_15_scaled)

# snow <- raster(file.path(GLs_folder, reg, "Snow_newValues_WGS84.tif"))
scaleAridity <- raster(file.path(GLs_folder, reg, "Aridity_FGRichnessScaled.tif"))
scaleAridity <- as.vector(scaleAridity)

ScalePETSD <- raster(file.path(GLs_folder, reg, "PETSD_FGRichnessScaled.tif"))
ScalePETSD <- as.vector(ScalePETSD)

scaleCLYPPT <- raster(file.path(GLs_folder, reg, "CLYPPT_FGRichnessCutScaled.tif"))
scaleCLYPPT <- as.vector(scaleCLYPPT)

scalePH <- raster(file.path(GLs_folder, reg, "PHIHOX_FGRichnessCutScaled.tif"))
scalePH <- as.vector(scalePH)

scaleSLTPPT <- raster(file.path(GLs_folder, reg, "SLTPPT_FGRichnessCutScaled.tif"))
scaleSLTPPT <- as.vector(scaleSLTPPT)




ESA <- raster(file.path(GLs_folder, reg, "ESA_newValuesCropped.tif"))
ESA <- as.vector(ESA)
keep <- c(60, 50, 70, 90, 110, 120, 130, 10, 12, 30)
ESA <- ifelse(ESA %in% keep, ESA, NA)
ESA <- as.factor(ESA)

newdat <- data.frame(ESA = ESA,
                     scaleSLTPPT = scaleSLTPPT,
                     scaleCLYPPT = scaleCLYPPT,
                     scalePH = scalePH,
                     ScalePETSD = ScalePETSD,
                     scaleAridity = scaleAridity, 
                     bio10_15_scaled = bio10_15_scaled,
                     bio10_1_scaled = bio10_1_scaled)


rm(list=c("bio10_1_scaled", "bio10_15_scaled",  "scaleAridity", "ScalePETSD",
          "scalePH", "scaleCLYPPT", "scaleSLTPPT","ESA"))

#############################################################

print("Predicting values...")
x <- split(newdat, (0:nrow(newdat) %/% 3000))  # modulo division

res <- c()

for(l in 1:length(x)){
  
  res <- c(res, predict(mod, x[[l]], re.form = NA))
  
}


length(res) == nrow(newdat)   

# need number of rows of the original raster
# The resolution
# the extent
# the coord.ref
# dimensions
# resol


print("Converting to raster...")
r <- matrix(res, nrow = dimensions[1], ncol = dimensions[2], byrow = TRUE)
r <- raster(r)

print("Adding in the raster information")

extent(r) <- exten
# ... and assign a projection
projection(r) <- coordred


# Save raster
print("Saving raster...")
r <- writeRaster(r,  filename=filename = file.path(savefolder, reg, "FGRichnessFinalRaster.tif"), format="GTiff", overwrite=TRUE)


print("Done!") 

