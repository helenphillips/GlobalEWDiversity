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
library(DHARMa)
source("Functions/FormatData.R")
source("Functions/lme4_ModellingFunctions.R")
source("Functions/ModelSimplification.R")
source("MEE3_1_sm_Appendix_S1/HighstatLib.R")
source("Functions/HypothesisTesting.R")

#################################################
# 2. Variables
#################################################


data_in <- "4_Data"
models <- "Models"

#################################################
# 3. Load data and models
#################################################
## Richness 
files <- list.files(file.path(data_in))
files <- files[grep("sitesRichness_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
richness <- read.csv(file.path(data_in, loadin))

# load(file.path(models, "richnessmodel_full.rds"))

## Biomass
files <- list.files(file.path(data_in))
files <- files[grep("sitesBiomass_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
biomass <- read.csv(file.path(data_in, loadin))

load(file.path(models, "biomassmodel_full.rds"))

## Abundance
files <- list.files(file.path(data_in))
files <- files[grep("sitesAbundance_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
abundance <- read.csv(file.path(data_in, loadin))

load(file.path(models, "abundancemodel_full.rds"))

#################################################
# 4. Abundance
#################################################

summary(abundance_model)


ESA <- "ESA"
mainEffects <- c("scalePH" , "scaleCLYPPT" ,"scaleSLTPPT" , "scaleCECSOL" , 
                   "scaleORCDRC" , "bio10_1_scaled" , "bio10_15_scaled" , "SnowMonths_cat" ,
                   "scaleAridity" , "ScalePETSD")
Climate <- c("bio10_1_scaled:bio10_15_scaled" , "bio10_1_scaled:SnowMonths_cat" ,  
  "bio10_1_scaled:scaleAridity" , "bio10_1_scaled:ScalePETSD" ,  
  "bio10_15_scaled:SnowMonths_cat" , "SnowMonths_cat:ScalePETSD")
  
Soil <- c("scalePH:scaleCLYPPT", "scalePH:scaleCECSOL", "scaleCLYPPT:scaleCECSOL", 
  "scaleSLTPPT:scaleCECSOL", "scaleCECSOL:scaleORCDRC")

WaterRetention <- c("scaleCLYPPT:bio10_15_scaled", "scaleCLYPPT:ScalePETSD")


HypothesisTesting(model = abundance_model, data = abundance, TestingGroups = c(ESA, Climate, Soil, WaterRetention))
