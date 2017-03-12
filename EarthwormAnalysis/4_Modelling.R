########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
library(lme4)
#################################################
# 2. Loading in variables
#################################################

data_in <-"3_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
loadin <- files[grep(date, files)]

rm(files)
rm(date)

if(!dir.exists("Models")){
  dir.create("Models")
}
models <- "Models"
#################################################
# 3. Load in data
#################################################

sites <- read.csv(file.path(data_in, loadin))
rm(loadin)


hist(sites$ph_new)

#################################################
# 4. Species Richness
#################################################

sites$scalePH <-scale(sites$ph_new)

sites_habitat <- droplevels(sites[sites$HabitatCover != "Unknown/Other",])

## There are some habitat covers that are not suitable for modelling at this stage
# sites_habitat <- droplevels(sites_habitat[!(sites_habitat$HabitatCover %in% c("Cropland/Other vegetation mosaic", "Paddy field", "Wetland")),])

sp_habitat <- glmer(NumberofSpecies ~ scalePH * HabitatCover + (1|Study_Name), data = sites_habitat, family = poisson)
summary(sp_habitat)
sp_habitat2 <- update(sp_habitat, .~. -scalePH:HabitatCover)
anova(sp_habitat, sp_habitat2) ## Significant at 0.01 level
plot(sp_habitat)
save(sp_habitat, file = file.path(models, "sp_habitat.rds"))



sites_lu <- droplevels(sites[sites$LandUse != "Unknown",])
sp_lu <- glmer(NumberofSpecies ~ scalePH * LandUse + (1|Study_Name), data = sites_lu, family = poisson)
summary(sp_lu)
sp_lu2 <- update(sp_lu, .~. -scalePH:LandUse)
anova(sp_lu, sp_lu2) ## Highly significant
save(sp_lu, file = file.path(models, "sp_landuse.rds"))
