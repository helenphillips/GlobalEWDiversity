##############################################
## -1. Editted colour picking function
##############################################

ColourPicker <- function(variable){
  if(any(levels(variable) %in% c("Pasture", "Primary vegetation", "Production - Arable"))){
    colos <- c(rep(NA, length(levels(variable))))
    
    tp <- data.frame(habitat = c("Pasture","Primary vegetation", "Production - Arable","Production - Crop plantations",
                                 "Production - Wood plantation", "Secondary vegetation", "Urban", "Unknown"),
                     colour =c("9BCD9B","1B5E20", "FFEB3B", "FF9800", "795548", "4CAF50",  "E65100",  "607D8B"))
  }
  colos <- tp$colour[match(levels(variable), tp$habitat)]
  return(colos)
}


########################################################
# 0. Set Working Directory
########################################################

setwd("C:/Users/hp39wasi/sWorm/UrbanforErin")



#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
library(lme4)
library(car)
library(Hmisc)

source("../EarthwormAnalysis/Functions/FormatData.R")
source("../EarthwormAnalysis/Functions/lme4_ModellingFunctions.R")
source("../EarthwormAnalysis/Functions/ModelSimplification.R")
source("../EarthwormAnalysis/Functions/Plots.R")
# source("../EarthwormAnalysis/Functions/ColourPicker.R")

#################################################
# 2. Loading in variables
#################################################

data_in <-"../EarthwormAnalysis/3_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

rm(files)
rm(date)

#################################################
# 2.5 Create folders
#################################################

if(!dir.exists("Data")){
  dir.create("Data")
}

data_out <- "Data"

if(!dir.exists("Models")){
  dir.create("Models")
}
models <- "Models"

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"
#################################################
# 3. Load in data
#################################################

sites <- read.csv(file.path(data_in, loadin))
# sites <- read.csv("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\3_Data\\Sites_2017-11-09.csv")
rm(loadin)
sites <- SiteLevels(sites) ## relevels all land use/habitat variables


##################################################
##
##################################################

table(sites$LandUse)
## 11 sites in urban
table(sites$ESA)

sites$Study_Name[sites$LandUse == "Urban"]
#################################################
## SPECIES RICHNESS MODELS
################################################
richness <- sites[complete.cases(sites$SpeciesRichness),] #3119
richness <- droplevels(richness[richness$LandUse != "Unknown",])
richness <- droplevels(richness[-which(richness$SpeciesRichness != round(richness$SpeciesRichness)),])

# richness <- richness[complete.cases(richness$scalePH),]

richness <- droplevels(richness[!(is.na(richness$PHIHOX)),])
# richness <- droplevels(richness[!(is.na(richness$bio10_2)),]) ## 2783

summary(richness$PHIHOX)
richness$scalePH <- as.vector(scale(richness$PHIHOX))
richness$scaleCLYPPT <- scale(richness$CLYPPT)
richness$scaleSLTPPT <- scale(richness$SLTPPT)
richness$scaleCECSOL <- scale(richness$CECSOL)
richness$scaleORCDRC <- scale(richness$ORCDRC)


# Three most important climate variables
richness$bio10_1_scaled <- scale(richness$bio10_1)
richness$bio10_4_scaled <- scale(richness$bio10_4)
# richness$bio10_7_scaled <- scale(richness$bio10_7)
richness$bio10_12_scaled <- scale(richness$bio10_12)
richness$bio10_15_scaled <- scale(richness$bio10_15)
## Save the data
write.csv(richness, file = file.path(data_out, paste("sitesRichness_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

r1 <- glmer(SpeciesRichness ~  LandUse + (scalePH  + 
                                        scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
            #  (bio10_1_scaled + bio10_4_scaled + 
             #    bio10_12_scaled + bio10_15_scaled)^2 + 
              #  SNDPPT # Not included, as the other two dictate the third
              
              #  (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
              
              # HabitatCover + 
              #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
              # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
              (1|file/Study_Name), data = richness, family = poisson,
            control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))


summary(r1)


richness_model <- modelSimplification(model = r1, data = richness, alpha = 0.05, optimizer = "bobyqa", Iters = 2e5)
save(richness_model, file = file.path(models, "richnessmodel_full.rds"))
# load(file.path(models, "richnessmodel_full.rds"))
richnessCols <- droplevels(ColourPicker(richness$LandUse))
jpeg(file = file.path(figures, "richness_landuse.jpg"), quality = 100, res = 200, width = 2000, height = 1500)
plotSingle(model= richness_model, 
           modelFixedEffs = c("LandUse", "scaleSLTPPT", "scaleCECSOL", "scaleORCDRC"),
           Effect1 = "LandUse", pt_cex = 2, pt_lwd = 3, axis_size = 1.2,
           responseVar = "SpeciesRichness", seMultiplier = 1, data = richness, cols = richnessCols, 
           legend.position, ylabel = "log-Species Richness", xlabel = "", otherContEffectsFun = "median")
dev.off()

###############################################################3
## TOTAL ABUNDANCE MODELS
###############################################

abundance <- sites[complete.cases(sites$Site_Abundance),]
abundance <- droplevels(abundance[abundance$LandUse != "Unknown",]) #3455

abundance <- droplevels(abundance[!(is.na(abundance$PHIHOX)),])

abundance$scalePH <- as.vector(scale(abundance$PHIHOX))
abundance$scaleCLYPPT <- scale(abundance$CLYPPT)
abundance$scaleSLTPPT <- scale(abundance$SLTPPT)
abundance$scaleCECSOL <- scale(abundance$CECSOL)
abundance$scaleORCDRC <- scale(abundance$ORCDRC)

## Save the data
write.csv(abundance, file = file.path(data_out, paste("sitesAbundance_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)


a1 <- lmer(logAbundance ~  LandUse + (scalePH  + 
                                    scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
            # (bio10_1_scaled + bio10_4_scaled +
            #   bio10_12_scaled + bio10_15_scaled)^2 +
             #  SNDPPT # Not included, as the other two dictate the third
             
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = abundance,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

plot(a1)
abundance_model <- modelSimplification(model = a1, data = abundance, alpha = 0.05, optimizer = "bobyqa", Iters = 2e5)
save(abundance_model, file = file.path(models, "abundancemodel_full.rds"))

jpeg(file = file.path(figures, "abundance_landuse.jpg"), quality = 100, res = 200, width = 2000, height = 1500)
abundanceCols <- droplevels(ColourPicker(abundance$LandUse))
plotSingle(model= abundance_model, 
           modelFixedEffs = c("LandUse", "scaleSLTPPT", "scaleCECSOL", "scaleORCDRC", "scalePH", "scaleCLYPPT"),
           Effect1 = "LandUse", 
           responseVar = "logAbundance", seMultiplier = 1, data = abundance, cols = abundanceCols, 
           legend.position, ylabel = "log-Abundance", xlabel = "", otherContEffectsFun = "median")
dev.off()



###########################################################3
## Sampling
############################################################

sites$samples <- paste(sites$Sampled_Area, sites$Sampled_Area_Unit)
table(sites$samples, sites$ESA)
