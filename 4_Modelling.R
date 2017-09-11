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


#plot(sp_habitat)
#save(sp_habitat, file = file.path(models, "sp_habitat.rds"))

#################################################
# 5. Biomass
#################################################
biomass <- sites[complete.cases(sites$Site_WetBiomass),]
biomass <- droplevels(biomass[biomass$LU_Mgmt != "Unknown",])
b1 <- lmer(log(Site_WetBiomass +1) ~ LU_Mgmt + scalePH + 
             scalePH:LU_Mgmt +
             # HabitatCover + 
            #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
            # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
            (1|file/Study_Name), data = biomass)

#################################################
# 6. Abundance
#################################################
hist(sites$Site_Abundance)
hist(log(sites$Site_Abundance + 1))
abundance <- sites[complete.cases(sites$Site_Abundance),]
abundance <- droplevels(abundance[abundance$LU_Mgmt != "Unknown",])

abundance$logAbundance <- log(abundance$Site_Abundance +1)

a1 <- lmer(logAbundance ~ LU_Mgmt + scalePH + 
             scalePH:LU_Mgmt + # HabitatCover + 
              # Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = abundance)

plot(a1)
#save(sp_habitat, file = file.path(models, "sp_habitat.rds"))


