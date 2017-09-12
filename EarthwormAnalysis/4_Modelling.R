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
library(car)
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

#################################################
# 4. Set reference levels
#################################################
sites$Management_System <- factor(sites$Management_System, levels = c("None",  "Annual crop", "Perennial crops", "Integrated systems",                   
                                                                      "Tree plantations", "Pastures (grazed lands)","Unknown"))

sites$LandUse <- factor(sites$LandUse, levels = c("Primary vegetation", "Secondary vegetation", "Pasture" ,
                                                  "Production - Arable", "Production - Crop plantations", 
                                                  "Production - Wood plantation",
                                                  "Urban", "Unknown"))

sites$HabitatCover <- factor(sites$HabitatCover, levels = c("Broadleaf deciduous forest", "Broadleaf evergreen forest",
                                                            "Needleleaf deciduous forest","Needleleaf evergreen forest",
                                                            "Mixed forest", "Tree open",
                                                            "Cropland","Cropland/Other vegetation mosaic",
                                                            "Herbaceous", "Herbaceous with spare tree/shrub",
                                                            "Shrub", "Sparse vegetation",
                                                            "Urban","Bare area (consolidated",
                                                            "Paddy field","Wetland", "Water bodies", "Unknown"))

sites$LU_Mgmt <- factor(sites$LU_Mgmt, levels = c( "Primary vegetation", "Secondary vegetation", "Annual crop", "Perennial crops",
                                                   "Integrated systems", "Tree plantations", "Pastures (grazed lands)", 
                                                   "Unknown", "Urban" ))

sites$intensity <- as.factor(sites$intensity)
#################################################
# 4. Species Richness
#################################################




#plot(sp_habitat)
#save(sp_habitat, file = file.path(models, "sp_habitat.rds"))

#################################################
# 5. Biomass
#################################################
biomass <- sites[complete.cases(sites$logBiomass),]
biomass <- droplevels(biomass[biomass$LU_Mgmt != "Unknown",])
b1 <- lmer(logBiomass ~ LU_Mgmt + scalePH + intensity +
             scalePH:LU_Mgmt +
             LU_Mgmt:intensity +
             # HabitatCover + 
            #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
            # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
            (1|file/Study_Name), data = biomass)

b2a <- update(b1, .~. -LU_Mgmt:intensity)
anova(b1, b2a) ## Not significant
b2b <- update(b1, .~. -LU_Mgmt:scalePH)
anova(b1, b2b) ##Significant

summary(b2a)
b3a <- update(b2a, .~. -LU_Mgmt:scalePH)
anova(b2a, b3a) ## Significant

b4a <- update(b2a, .~. -intensity)
anova(b2a, b4a) ## Significant

######
## b2a
######


#################################################
# 6. Abundance
#################################################
hist(sites$Site_Abundance)
hist(sites$logAbundance)
abundance <- sites[complete.cases(sites$Site_Abundance),]
abundance <- droplevels(abundance[abundance$LU_Mgmt != "Unknown",])

a1 <- lmer(logAbundance ~ LU_Mgmt + scalePH + intensity +
             scalePH:LU_Mgmt + LU_Mgmt:intensity + # HabitatCover + 
              # Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = abundance)

plot(a1)
#save(sp_habitat, file = file.path(models, "sp_habitat.rds"))


a2a <- update(a1, .~. -LU_Mgmt:intensity)
anova(a1, a2a) ## Not quite significant
a2b <- update(a1, .~. -LU_Mgmt:scalePH)
anova(a1, a2b) # Significant

a3a <- update(a2a, .~. -LU_Mgmt:scalePH)
anova(a2a, a3a) ## Significant

a4a <- update(a2a, .~. -intensity)
anova(a2a, a4a) ## Significant

####
## a2a
####