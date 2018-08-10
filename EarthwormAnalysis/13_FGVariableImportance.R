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

load(file.path(models, "richnessmodel_epifunctionalgroups.rds"))
load(file.path(models, "richnessmodel_endofunctionalgroups.rds"))
load(file.path(models, "richnessmodel_anefunctionalgroups.rds"))

## Abundance 
files <- list.files(file.path(data_in))
files <- files[grep("sitesFGAbundance_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
abundance <- read.csv(file.path(data_in, loadin))

load(file.path(models, "abundancemodel_epifunctionalgroups.rds"))
load(file.path(models, "abundancemodel_endofunctionalgroups.rds"))
load(file.path(models, "abundancemodel_anefunctionalgroups.rds"))

## Biomass 
files <- list.files(file.path(data_in))
files <- files[grep("sitesFGBiomass_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
biomass <- read.csv(file.path(data_in, loadin))

load(file.path(models, "biomassmodel_epifunctionalgroups.rds"))
load(file.path(models, "biomassmodel_endofunctionalgroups.rds"))
load(file.path(models, "biomassmodel_anefunctionalgroups.rds"))

###############################################################
##  BIOMASS - with split groups
###############################################################

# BIOMASS
# epi
ESA <- "ESA"
epi_Temperature <- c("bio10_4_scaled", "ScalePETSD")
epi_Precip <- c("bio10_12_scaled", "bio10_15_scaled", "SnowMonths_cat")
epi_Soil <- c("scalePH","scaleORCDRC","scaleCLYPPT", "scaleSLTPPT")
epi_WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_12_scaled", "bio10_15_scaled", "ScalePETSD")

groups <- list(
  ESA = ESA,
  epi_Temperature = epi_Temperature,
  epi_Precip = epi_Precip,
  epi_Soil = epi_Soil,
  epi_WaterRetention = epi_WaterRetention
)


epi_biomass <- biomass[grep("Epi", biomass$variable),]
epi_biomass <- scaleVariables(epi_biomass)
epi_biomass$value[which(epi_biomass$value < 0)] <- 0
epi_biomass$logValue <- log(epi_biomass$value + 1)

epi_biomass_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "bio10_4_scaled", "bio10_12_scaled",  
                         "bio10_15_scaled", "SnowMonths_cat", "ScalePETSD", "ESA")


epi_bioM_rf <- randomForest(y = epi_biomass$logValue, x = epi_biomass[,names(epi_biomass) %in% epi_biomass_mainEffects], 
                        ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(epi_bioM_rf, type=1)
varImpPlot(epi_bioM_rf, type=2)

epi_biomass_import_split <- group.importance(epi_bioM_rf, groups = groups) 

# endo
endo_Temperature <- c("bio10_4_scaled", "ScalePETSD")
endo_Precip <- c("bio10_12_scaled", "bio10_15_scaled", "SnowMonths_cat")
endo_Soil <- c("scalePH","scaleORCDRC","scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL")
endo_WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_12_scaled", "bio10_15_scaled", "ScalePETSD")

groups <- list(
  ESA = ESA,
  endo_Temperature = endo_Temperature,
  endo_Precip = endo_Precip,
  endo_Soil = endo_Soil,
  endo_WaterRetention = endo_WaterRetention
)


endo_biomass <- biomass[grep("Endo", biomass$variable),]
endo_biomass <- scaleVariables(endo_biomass)
endo_biomass$value[which(endo_biomass$value < 0)] <- 0
endo_biomass$logValue <- log(endo_biomass$value + 1)

endo_biomass_mainEffects <- c("scalePH", "scaleCECSOL", "scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "bio10_4_scaled", "bio10_12_scaled",  
                             "bio10_15_scaled", "SnowMonths_cat", "ScalePETSD", "ESA")


endo_bioM_rf <- randomForest(y = endo_biomass$logValue, x = endo_biomass[,names(endo_biomass) %in% endo_biomass_mainEffects], 
                            ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(endo_bioM_rf, type=1)
varImpPlot(endo_bioM_rf, type=2)

endo_biomass_import_split <- group.importance(endo_bioM_rf, groups = groups) 

############
# Biomass figure
############
## 5 is most important, 1 is least important

epi_biomass_import_split
epi_biomass_order <- c(1, 5, 4, 2, 3)
endo_biomass_import_split
endo_biomass_order <- c(5,4,2,1,3)

ane_biomass_order <- c(NA, NA, NA, NA, NA)

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
ESA <- "ESA"
epi_Temperature <- c("bio10_4_scaled", "ScalePET")
epi_Precip <- c("bio10_15_scaled", "SnowMonths_cat", "scaleAridity")
epi_Soil <- c("scalePH","scaleORCDRC","scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL")
epi_WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_15_scaled", "ScalePET", "scaleAridity")

groups <- list(
  ESA = ESA,
  epi_Temperature = epi_Temperature,
  epi_Precip = epi_Precip,
  epi_Soil = epi_Soil,
  epi_WaterRetention = epi_WaterRetention
)


epi_abundance <- abundance[grep("Epi", abundance$variable),]
epi_abundance <- scaleVariables(epi_abundance)
epi_abundance$logValue <- log(epi_abundance$value + 1)

epi_abundance_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "scaleCECSOL", "bio10_4_scaled", 
                             "bio10_15_scaled", "SnowMonths_cat", "scaleAridity", "ScalePET", "ESA")


epi_abund_rf <- randomForest(y = epi_abundance$logValue, x = epi_abundance[,names(epi_abundance) %in% epi_abundance_mainEffects], 
                            ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(epi_abund_rf, type=1)
varImpPlot(epi_abund_rf, type=2)

epi_abundance_import_split <- group.importance(epi_abund_rf, groups = groups) 

# endo
endo_Temperature <- c("bio10_4_scaled", "ScalePET")
endo_Precip <- c("bio10_15_scaled", "SnowMonths_cat", "scaleAridity")
endo_Soil <- c("scalePH","scaleORCDRC","scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL")
endo_WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "scaleAridity", "bio10_15_scaled", "ScalePET")

groups <- list(
  ESA = ESA,
  endo_Temperature = endo_Temperature,
  endo_Precip = endo_Precip,
  endo_Soil = endo_Soil,
  endo_WaterRetention = endo_WaterRetention
)


endo_abundance <- abundance[grep("Endo", abundance$variable),]
endo_abundance <- scaleVariables(endo_abundance)
endo_abundance$value[which(endo_abundance$value < 0)] <- 0
endo_abundance$logValue <- log(endo_abundance$value + 1)

endo_abundance_mainEffects <- c("scalePH", "scaleCECSOL", "scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "bio10_4_scaled", "bio10_15_scaled",  
                              "SnowMonths_cat", "ScalePET", "scaleAridity", "ESA")


endo_abund_rf <- randomForest(y = endo_abundance$logValue, x = endo_abundance[,names(endo_abundance) %in% endo_abundance_mainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(endo_abund_rf, type=1)
varImpPlot(endo_abund_rf, type=2)

endo_abundance_import_split <- group.importance(endo_abund_rf, groups = groups) 

# ane
ane_Temperature <- c("bio10_4_scaled", "ScalePET")
ane_Precip <- c("bio10_15_scaled", "SnowMonths_cat", "scaleAridity")
ane_Soil <- c("scalePH","scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL")
ane_WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "scaleAridity", "bio10_15_scaled", "ScalePET")

groups <- list(
  ESA = ESA,
  ane_Temperature = ane_Temperature,
  ane_Precip = ane_Precip,
  ane_Soil = ane_Soil,
  ane_WaterRetention = ane_WaterRetention
)


ane_abundance <- abundance[grep("Ane", abundance$variable),]
ane_abundance <- scaleVariables(ane_abundance)
ane_abundance$value[which(ane_abundance$value < 0)] <- 0
ane_abundance$logValue <- log(ane_abundance$value + 1)

ane_abundance_mainEffects <- c("scalePH", "scaleCECSOL", "scaleCLYPPT", "scaleSLTPPT", "bio10_4_scaled", "bio10_15_scaled",  
                                "SnowMonths_cat", "ScalePET", "scaleAridity", "ESA")


ane_abund_rf <- randomForest(y = ane_abundance$logValue, x = ane_abundance[,names(ane_abundance) %in% ane_abundance_mainEffects], 
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
ane_abundance_order <- c(1,5,2,4,3)

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
# epi
ESA <- "ESA"
epi_Temperature <- c("ScalePET", "ScalePETSD")
epi_Precip <- c("bio10_12_scaled", "SnowMonths_cat")
epi_Soil <- c("scalePH","scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL")
epi_WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_12_scaled", "ScalePET", "ScalePETSD")

groups <- list(
  ESA = ESA,
  epi_Temperature = epi_Temperature,
  epi_Precip = epi_Precip,
  epi_Soil = epi_Soil,
  epi_WaterRetention = epi_WaterRetention
)


epi_richness <- richness[grep("Epi", richness$variable),]
epi_richness <- scaleVariables(epi_richness)

epi_richness_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL",  
                               "bio10_12_scaled", "SnowMonths_cat", "ScalePET","ScalePETSD", "ESA")


epi_richness_rf <- randomForest(y = epi_richness$value, x = epi_richness[,names(epi_richness) %in% epi_richness_mainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(epi_richness_rf, type=1)
varImpPlot(epi_richness_rf, type=2)

epi_richness_import_split <- group.importance(epi_richness_rf, groups = groups) 

# endo
endo_Temperature <- c("bio10_4_scaled")
endo_Precip <- c("bio10_15_scaled")
endo_Soil <- c("scalePH", "scaleSLTPPT", "scaleCECSOL")
endo_WaterRetention <- c("scaleSLTPPT", "bio10_15_scaled")

groups <- list(
  # ESA = ESA,
  endo_Temperature = endo_Temperature,
  endo_Precip = endo_Precip,
  endo_Soil = endo_Soil,
  endo_WaterRetention = endo_WaterRetention
)


endo_richness <- richness[grep("Endo", richness$variable),]
endo_richness <- scaleVariables(endo_richness)

endo_richness_mainEffects <- c("scalePH", "scaleCECSOL", "scaleSLTPPT",  "bio10_4_scaled", "bio10_15_scaled")


endo_richness_rf <- randomForest(y = endo_richness$value, x = endo_richness[,names(endo_richness) %in% endo_richness_mainEffects], 
                              ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(endo_richness_rf, type=1)
varImpPlot(endo_richness_rf, type=2)

endo_richness_import_split <- group.importance(endo_richness_rf, groups = groups) 

# ane
ane_Temperature <- c("bio10_1_scaled")
ane_Precip <- c("bio10_12_scaled", "SnowMonths_cat")
ane_Soil <- c("scalePH","scaleSLTPPT")
ane_WaterRetention <- c("scaleSLTPPT", "bio10_12_scaled")

groups <- list(
  # ESA = ESA,
  ane_Temperature = ane_Temperature,
  ane_Precip = ane_Precip,
  ane_Soil = ane_Soil,
  ane_WaterRetention = ane_WaterRetention
)


ane_richness <- richness[grep("Ane", richness$variable),]
ane_richness <- scaleVariables(ane_richness)

ane_richness_mainEffects <- c("scalePH", "scaleSLTPPT", "bio10_1_scaled", "bio10_12_scaled",  
                               "SnowMonths_cat")


ane_richness_rf <- randomForest(y = ane_richness$value, x = ane_richness[,names(ane_richness) %in% ane_richness_mainEffects], 
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
epi_richness_order <- c(1, 5, 2, 4, 3)
endo_richness_import_split
endo_richness_order <- c(NA, 5,2,4,3)
ane_richness_import_split
ane_richness_order <- c(NA,5,3,2,4)

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


