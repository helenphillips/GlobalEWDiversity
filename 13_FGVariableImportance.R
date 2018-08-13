########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################

library(lme4)
library(randomForest)
library(ggplot2)
library(reshape)
library(viridis)
library(RColorBrewer)

source("Functions/FormatData.R")
source("Functions/lme4_ModellingFunctions.R")
source("Functions/ModelSimplification.R")
source("MEE3_1_sm_Appendix_S1/HighstatLib.R")
source("Functions/RF_VariableSetImportance.R")
source("Functions/Plots.R")

#################################################
# 2. Variables
#################################################
data_in <- "11_Data"
models <- "Models"
figures <- "Figures"

#################################################
# 3. Load data and models
#################################################
## Richness 
files <- list.files(file.path(data_in))
files <- files[grep("sitesFGRichness_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
richness <- read.csv(file.path(data_in, loadin))


## Abundance 
files <- list.files(file.path(data_in))
files <- files[grep("sitesFGAbundance_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
abundance <- read.csv(file.path(data_in, loadin))


## Biomass 
files <- list.files(file.path(data_in))
files <- files[grep("sitesFGBiomass_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
biomass <- read.csv(file.path(data_in, loadin))

############################################################
# Load in main models
###########################################################
# Using the variables from these
load(file.path(models, "richnessmodel_full.rds"))
load(file.path(models, "biomassmodel_full.rds"))
load(file.path(models, "abundancemodel_full.rds"))


###############################################################
##  BIOMASS - with split groups
###############################################################
# logBiomass ~ scalePH + scaleCLYPPT + scaleSLTPPT + scaleORCDRC +  
# scaleCECSOL + bio10_12_scaled + bio10_15_scaled + scalePH:scaleCLYPPT +  
#  scalePH:scaleSLTPPT + scalePH:scaleORCDRC + scalePH:scaleCECSOL +  
#  scaleCLYPPT:scaleORCDRC + scaleCLYPPT:scaleCECSOL + scaleSLTPPT:scaleCECSOL +  
#  scaleORCDRC:scaleCECSOL + scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +  
#  scaleCLYPPT:bio10_15_scaled + ESA + SnowMonths_cat


# BIOMASS
biomass_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "scaleCECSOL",
                         "bio10_12_scaled", "bio10_15_scaled", "SnowMonths_cat", "ESA")

ESA <- "ESA"
# biomass_Temperature <- c("bio10_4_scaled", "ScalePETSD")
biomass_Precip <- c("bio10_12_scaled", "bio10_15_scaled", "SnowMonths_cat")
biomass_Soil <- c("scalePH","scaleORCDRC","scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL")
biomass_WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_12_scaled", "bio10_15_scaled")

groups <- list(
  ESA = ESA,
  # biomass_Temperature = biomass_Temperature,
  biomass_Precip = biomass_Precip,
  biomass_Soil = biomass_Soil,
  biomass_WaterRetention = biomass_WaterRetention
)


epi_biomass <- biomass[grep("Epi", biomass$variable),]
epi_biomass <- scaleVariables(epi_biomass)
epi_biomass$value[which(epi_biomass$value < 0)] <- 0
epi_biomass$logValue <- log(epi_biomass$value + 1)

epi_bioM_rf <- randomForest(y = epi_biomass$logValue, x = epi_biomass[,names(epi_biomass) %in% biomass_mainEffects], 
                        ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(epi_bioM_rf, type=1)
varImpPlot(epi_bioM_rf, type=2)

epi_biomass_import_split <- group.importance(epi_bioM_rf, groups = groups) 

# endo
endo_biomass <- biomass[grep("Endo", biomass$variable),]
endo_biomass <- scaleVariables(endo_biomass)
endo_biomass$value[which(endo_biomass$value < 0)] <- 0
endo_biomass$logValue <- log(endo_biomass$value + 1)

endo_bioM_rf <- randomForest(y = endo_biomass$logValue, x = endo_biomass[,names(endo_biomass) %in% biomass_mainEffects], 
                            ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(endo_bioM_rf, type=1)
varImpPlot(endo_bioM_rf, type=2)

endo_biomass_import_split <- group.importance(endo_bioM_rf, groups = groups) 

## Anecics
ane_biomass <- biomass[grep("Ane", biomass$variable),]
ane_biomass <- scaleVariables(ane_biomass)
ane_biomass$value[which(ane_biomass$value < 0)] <- 0
ane_biomass$logValue <- log(ane_biomass$value + 1)

ane_bioM_rf <- randomForest(y = ane_biomass$logValue, x = ane_biomass[,names(ane_biomass) %in% biomass_mainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(ane_bioM_rf, type=1)
varImpPlot(ane_bioM_rf, type=2)

ane_biomass_import_split <- group.importance(ane_bioM_rf, groups = groups) 


############
# Biomass figure
############
## 5 is most important, 1 is least important
# NA is 1 
epi_biomass_import_split
epi_biomass_order <- c(2, NA, 5, 3, 4)
endo_biomass_import_split
endo_biomass_order <- c(5,NA,4,2,3)
ane_biomass_import_split
ane_biomass_order <- c(5, NA, 3, 4, 2)

a <- matrix(rep(NA, length = 5*3), nrow = 3, ncol = 5)
colnames(a) <- c("ESA", "Temperature", "Precipitation", "Soil", "Water Retention")
a[1,] <- epi_biomass_order
a[2,] <- endo_biomass_order
a[3,] <- ane_biomass_order
rownames(a) <- c("Epigeics", "Endogeics", "Anecics")

dat <- melt(a)

dat$X1 <- factor(dat$X1, levels = c("Anecics", "Endogeics","Epigeics"))
dat$X2 <- factor(dat$X2, levels = c("ESA","Soil","Precipitation", "Temperature","Water Retention"))

jpeg(file = file.path(figures, "variableImportance_BiomassFGs.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
p <- VariableImportancePlot(dat, lowColour = "#BCBDDC", highColour = "#25004b", yLab = "Functional Group Model")
p
dev.off()

###############################################################
##  ABUNDANCE - with split groups
###############################################################
# epi
# logAbundance ~ scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL +  
#  scaleORCDRC + bio10_1_scaled + bio10_15_scaled + SnowMonths_cat +  
#  scaleAridity + ScalePETSD + scalePH:scaleCLYPPT + scalePH:scaleCECSOL +  
#  scaleCLYPPT:scaleCECSOL + scaleSLTPPT:scaleCECSOL + scaleCECSOL:scaleORCDRC +  
#  bio10_1_scaled:bio10_15_scaled + bio10_1_scaled:SnowMonths_cat +  
#  bio10_1_scaled:scaleAridity + bio10_1_scaled:ScalePETSD +  
#  bio10_15_scaled:SnowMonths_cat + SnowMonths_cat:ScalePETSD +  
#  scaleCLYPPT:bio10_15_scaled + scaleCLYPPT:ScalePETSD + ESA +  

abundance_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "scaleCECSOL", 
                           "bio10_1_scaled", "bio10_15_scaled", "SnowMonths_cat", "scaleAridity", 
                           "ScalePETSD", "ESA")


ESA <- "ESA"
abundance_Temperature <- c("bio10_1_scaled", "ScalePETSD")
abundance_Precip <- c("bio10_15_scaled", "SnowMonths_cat", "scaleAridity")
abundance_Soil <- c("scalePH","scaleORCDRC","scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL")
abundance_WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_15_scaled", "ScalePETSD", "scaleAridity")

groups <- list(
  ESA = ESA,
  abundance_Temperature = abundance_Temperature,
  abundance_Precip = abundance_Precip,
  abundance_Soil = abundance_Soil,
  abundance_WaterRetention = abundance_WaterRetention
)


epi_abundance <- abundance[grep("Epi", abundance$variable),]
epi_abundance <- scaleVariables(epi_abundance)
epi_abundance$logValue <- log(epi_abundance$value + 1)

epi_abund_rf <- randomForest(y = epi_abundance$logValue, x = epi_abundance[,names(epi_abundance) %in% abundance_mainEffects], 
                            ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(epi_abund_rf, type=1)
varImpPlot(epi_abund_rf, type=2)

epi_abundance_import_split <- group.importance(epi_abund_rf, groups = groups) 

# endo
endo_abundance <- abundance[grep("Endo", abundance$variable),]
endo_abundance <- scaleVariables(endo_abundance)
endo_abundance$value[which(endo_abundance$value < 0)] <- 0
endo_abundance$logValue <- log(endo_abundance$value + 1)

endo_abund_rf <- randomForest(y = endo_abundance$logValue, x = endo_abundance[,names(endo_abundance) %in% abundance_mainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(endo_abund_rf, type=1)
varImpPlot(endo_abund_rf, type=2)

endo_abundance_import_split <- group.importance(endo_abund_rf, groups = groups) 

# ane
ane_abundance <- abundance[grep("Ane", abundance$variable),]
ane_abundance <- scaleVariables(ane_abundance)
ane_abundance$value[which(ane_abundance$value < 0)] <- 0
ane_abundance$logValue <- log(ane_abundance$value + 1)

ane_abund_rf <- randomForest(y = ane_abundance$logValue, x = ane_abundance[,names(ane_abundance) %in% abundance_mainEffects], 
                              ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(ane_abund_rf, type=1)
varImpPlot(ane_abund_rf, type=2)

ane_abundance_import_split <- group.importance(ane_abund_rf, groups = groups) 

############
# Abundance figure
############
## 5 is most important, 1 is least important

epi_abundance_import_split
epi_abundance_order <- c(1, 5, 2, 4, 3)
endo_abundance_import_split
endo_abundance_order <- c(5,4,1,3,2)
ane_abundance_import_split
ane_abundance_order <- c(1,5,4,3, 2)

a <- matrix(rep(NA, length = 5*3), nrow = 3, ncol = 5)
colnames(a) <- c("ESA", "Temperature", "Precipitation", "Soil", "Water Retention")
a[1,] <- epi_abundance_order
a[2,] <- endo_abundance_order
a[3,] <- ane_abundance_order
rownames(a) <- c("Epigeics", "Endogeics", "Anecics")

dat <- melt(a)

dat$X1 <- factor(dat$X1, levels = c("Anecics", "Endogeics","Epigeics"))
dat$X2 <- factor(dat$X2, levels = c("ESA","Soil","Precipitation", "Temperature","Water Retention"))

jpeg(file = file.path(figures, "variableImportance_AbundanceFGs.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
p <- VariableImportancePlot(dat, lowColour = "#BCBDDC", highColour = "#25004b", yLab = "Functional Group Model")
p
dev.off()

###############################################################
##  SPECIES RICHNESS - with split groups
###############################################################

# SpeciesRichness ~ scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL +  
#  scaleORCDRC + bio10_4_scaled + bio10_15_scaled + SnowMonths_cat +  
#  scaleAridity + ScalePET + scalePH:scaleCECSOL + scalePH:scaleORCDRC +  
#  scaleCLYPPT:scaleCECSOL + scaleSLTPPT:scaleORCDRC + scaleCECSOL:scaleORCDRC +  
#  bio10_4_scaled:bio10_15_scaled + bio10_4_scaled:SnowMonths_cat +  
#  bio10_15_scaled:SnowMonths_cat + bio10_15_scaled:scaleAridity +  
#  bio10_15_scaled:ScalePET + scaleCLYPPT:bio10_15_scaled +  
#  scaleSLTPPT:ScalePET + scaleSLTPPT:scaleAridity + ESA + 

# epi

ESA <- "ESA"
richness_Temperature <- c("bio10_4_scaled", "ScalePET")
richness_Precip <- c("bio10_15_scaled", "SnowMonths_cat", "scaleAridity")
richness_Soil <- c("scalePH","scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL", "scaleORCDRC")
richness_WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_15_scaled", "ScalePET", "scaleAridity")

groups <- list(
  ESA = ESA,
  richness_Temperature = richness_Temperature,
  richness_Precip = richness_Precip,
  richness_Soil = richness_Soil,
  richness_WaterRetention = richness_WaterRetention
)

richness_mainEffects <- c(ESA, richness_Temperature, richness_Precip, richness_Soil)


epi_richness <- richness[grep("Epi", richness$variable),]
epi_richness <- scaleVariables(epi_richness)

epi_richness_rf <- randomForest(y = epi_richness$value, x = epi_richness[,names(epi_richness) %in% richness_mainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(epi_richness_rf, type=1)
varImpPlot(epi_richness_rf, type=2)

epi_richness_import_split <- group.importance(epi_richness_rf, groups = groups) 

# endo
endo_richness <- richness[grep("Endo", richness$variable),]
endo_richness <- scaleVariables(endo_richness)

endo_richness_rf <- randomForest(y = endo_richness$value, x = endo_richness[,names(endo_richness) %in% richness_mainEffects], 
                              ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(endo_richness_rf, type=1)
varImpPlot(endo_richness_rf, type=2)

endo_richness_import_split <- group.importance(endo_richness_rf, groups = groups) 

# ane
ane_richness <- richness[grep("Ane", richness$variable),]
ane_richness <- scaleVariables(ane_richness)

ane_richness_rf <- randomForest(y = ane_richness$value, x = ane_richness[,names(ane_richness) %in% richness_mainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(ane_richness_rf, type=1)
varImpPlot(ane_richness_rf, type=2)

ane_richness_import_split <- group.importance(ane_richness_rf, groups = groups) 

############
# Species Richness figure
############
## 5 is most important, 1 is least important
# ESA, Temperature, Precip, Soil, WaterRetention

epi_richness_import_split
epi_richness_order <- c(1, 4, 2, 5, 3)
endo_richness_import_split
endo_richness_order <- c(1, 5,2,4,3)
ane_richness_import_split
ane_richness_order <- c(1,5,4,2,3)

a <- matrix(rep(NA, length = 5*3), nrow = 3, ncol = 5)
colnames(a) <- c("ESA", "Temperature", "Precipitation", "Soil", "Water Retention")
a[1,] <- epi_richness_order
a[2,] <- endo_richness_order
a[3,] <- ane_richness_order
rownames(a) <- c("Epigeics", "Endogeics", "Anecics")

dat <- melt(a)

dat$X1 <- factor(dat$X1, levels = c("Anecics", "Endogeics","Epigeics"))
dat$X2 <- factor(dat$X2, levels = c("ESA","Soil","Precipitation", "Temperature","Water Retention"))


jpeg(file = file.path(figures, "variableImportance_RichnessFGs.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
p <- VariableImportancePlot(dat, lowColour = "#BCBDDC", highColour = "#25004b", yLab = "Functional Group Model")
p
dev.off()


