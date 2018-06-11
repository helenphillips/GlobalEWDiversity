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
# source("Functions/HypothesisTesting.R")
source("Functions/RF_VariableSetImportance.R")

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

load(file.path(models, "richnessmodel_full.rds"))

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

Climate <- c("bio10_1_scaled", "bio10_15_scaled" , "SnowMonths_cat" ,  
  "scaleAridity" , "ScalePETSD")
  
Soil <- c("scalePH", "scaleCLYPPT", "scaleCECSOL", 
  "scaleSLTPPT", "scaleORCDRC")

WaterRetention <- c("scaleCLYPPT", "bio10_15_scaled", "ScalePETSD")

groups <- list(
  ESA = ESA,
  Climate = Climate,
  Soil = Soil,
  WaterRetention = WaterRetention
)

## Randomforests can usually still find the interaction
# https://stats.stackexchange.com/questions/157665/including-interaction-terms-in-random-forest?rq=1
# https://stats.stackexchange.com/questions/201893/how-to-include-an-interaction-term-in-a-random-forest-model
## So maybe still comparable to the final regression??

library(randomForest)

## Abundance

AbundmainEffects <- c("scalePH" , "scaleCLYPPT" ,"scaleSLTPPT" , "scaleCECSOL" , 
                 "scaleORCDRC" , "bio10_1_scaled" , "bio10_15_scaled" , "SnowMonths_cat" ,
                 "scaleAridity" , "ScalePETSD", "ESA")

abundance_rf <- randomForest(y = abundance$logAbundance, x = abundance[,names(abundance) %in% AbundmainEffects], 
                  ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(abundance_rf, type=1) # mean decrease in accuracy
# "The accuracy one tests to see how worse the model performs without each variable,
# so a high decrease in accuracy would be expected for very predictive variables. 
# The Gini one (type 2) digs into the mathematics behind decision trees, but essentially 
# measures how pure the nodes are at the end of the tree. Again it tests to see the 
# result if each variable is taken out and a high score means the variable was important."
# http://trevorstephens.com/post/73770963794/titanic-getting-started-with-r-part-5-random
varImpPlot(abundance_rf, type=2)

## Variable importance across a set
## USing GINI accuracy
abundance_import <- group.importance(abundance_rf, groups) 

############
# Richness


richness_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL",  
  "scaleORCDRC", "bio10_4_scaled", "bio10_15_scaled", "SnowMonths_cat",  
  "scaleAridity", "ScalePET", "ESA")


ESA <- "ESA"

Climate <- c("bio10_4_scaled", "bio10_15_scaled", "SnowMonths_cat",  
             "scaleAridity", "ScalePET")

Soil <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL", "scaleORCDRC")

WaterRetention <- c("scaleCLYPPT","scaleSLTPPT", "bio10_15_scaled", "scaleAridity", "ScalePET")

groups <- list(
  ESA = ESA,
  Climate = Climate,
  Soil = Soil,
  WaterRetention = WaterRetention
)


spR_rf <- randomForest(y = richness$SpeciesRichness, x = richness[,names(richness) %in% richness_mainEffects], 
                  ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(spR_rf, type=1)
richness_import <- group.importance(spR_rf, groups) 

############################
## Biomass


biomass_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "scaleCECSOL", "bio10_12_scaled",  
  "bio10_15_scaled", "SnowMonths_cat", "ESA")

ESA <- "ESA"

Climate <- c("bio10_12_scaled", "bio10_15_scaled", "SnowMonths_cat")

Soil <- c("scalePH","scaleORCDRC", "scaleCECSOL","scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC")

WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_12_scaled", "bio10_15_scaled")

groups <- list(
  ESA = ESA,
  Climate = Climate,
  Soil = Soil,
  WaterRetention = WaterRetention
)

bioM_rf <- randomForest(y = biomass$logBiomass, x = biomass[,names(biomass) %in% biomass_mainEffects], 
                       ntree=501, importance=TRUE, proximity = TRUE)
varImpPlot(bioM_rf, type=1)
varImpPlot(bioM_rf, type=2)
biomass_import <- group.importance(bioM_rf, groups) 

#######################################
## PLOTS
#######################################

a <- matrix(rep(NA, length = 4*3), nrow = 3, ncol = 4)
colnames(a) <- names(groups)
a[1,] <- as.vector(OrderImportance(richness_import))
a[2,] <- as.vector(OrderImportance(abundance_import))
a[3,] <- as.vector(OrderImportance(biomass_import))
## 4 is most important, 1 is least important
rownames(a) <- c("SpeciesRichness", "Abundance", "Biomass")




library(reshape)
dat <- melt(a)

library(ggplot2)
p <- ggplot(data =  dat, aes(x = X2, y = X1)) +
  geom_tile(aes(fill = value), colour = "white") +
  geom_text(aes(label = sprintf("%1.2f",value)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "grey28")
p
