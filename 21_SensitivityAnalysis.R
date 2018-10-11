if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}
#################################################
# 1. Libraries
#################################################
library(lme4)

createSplits <- function(dat, kfold = 10){
  rows <- nrow(dat)
  rows <- sample(rows, size = length(1:rows))
  
  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  
  splits <- chunk(rows, kfold)
  
  return(splits)
}
#################################################
# 2. Load data
#################################################

data_in <-"3.5_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

rm(files)
rm(date)

sites <- read.csv(file.path(data_in, loadin))

data_in <-"4_Data"
## Richness 
files <- list.files(file.path(data_in))
files <- files[grep("sitesRichness_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
richness <- read.csv(file.path(data_in, loadin))

## Biomass
files <- list.files(file.path(data_in))
files <- files[grep("sitesBiomass_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
biomass <- read.csv(file.path(data_in, loadin))

## Abundance
files <- list.files(file.path(data_in))
files <- files[grep("sitesAbundance_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
abundance <- read.csv(file.path(data_in, loadin))

#################################################
# 3. Create directories
#################################################

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"


#################################################
# 4. Load in models
#################################################

models <- "Models"


load(file.path(models, "richnessmodel.rds"))
load(file.path(models, "biomassmodel_full.rds"))
load(file.path(models, "abundancemodel_full.rds"))

#################################################
# 5. Species Richness
################################################
richnessData <- richness_model@frame


########
# K-Fold Cross validation
########
k_fold <- 10
splits <- createSplits(richnessData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  
  rows <- as.vector(unlist(splits[k]))
  testData <- richnessData[rows,]
  bankData <- richnessData[-rows,]
  
  mod <- glmer(SpeciesRichness ~ scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL +  
                scaleORCDRC + bio10_7_scaled + bio10_15_scaled + SnowMonths_cat +  
                scaleAridity + ScalePET + scalePH:scaleCECSOL + scalePH:scaleORCDRC +  
                scaleCLYPPT:scaleCECSOL + scaleSLTPPT:scaleORCDRC + scaleCECSOL:scaleORCDRC +  
                bio10_7_scaled:bio10_15_scaled + bio10_7_scaled:SnowMonths_cat +  
                bio10_7_scaled:scaleAridity + bio10_15_scaled:SnowMonths_cat +  
                bio10_15_scaled:scaleAridity + bio10_15_scaled:ScalePET +  
                SnowMonths_cat:scaleAridity + scaleCLYPPT:bio10_15_scaled +  
                scaleSLTPPT:ScalePET + scaleSLTPPT:scaleAridity + ESA + scaleElevation + 
                (1 | file/Study_Name), data = bankData, family = poisson,
              control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))
  
  testData$Predicted <- (predict(mod, testData,  re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$SpeciesRichness, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)
plot(df$predicted ~ df$observed)
abline(0, 1) 


#################################################
# 5. Abundance
################################################

abundanceData <- abundance_model@frame


########
# K-Fold Cross validation
########
  k_fold <- 10
splits <- createSplits(abundanceData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  
  rows <- as.vector(unlist(splits[k]))
  testData <- abundanceData[rows,]
  bankData <- abundanceData[-rows,]
  
  mod <- lmer(logAbundance ~ scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL +  
                scaleORCDRC + bio10_7_scaled + bio10_15_scaled + SnowMonths_cat +  
                scaleAridity + ScalePET + scalePH:scaleSLTPPT + scalePH:scaleCECSOL +  
                scalePH:scaleORCDRC + scaleCLYPPT:scaleCECSOL + scaleCLYPPT:scaleORCDRC +  
                scaleCECSOL:scaleORCDRC + bio10_7_scaled:bio10_15_scaled +  
                bio10_7_scaled:SnowMonths_cat + bio10_7_scaled:scaleAridity +  
                bio10_15_scaled:SnowMonths_cat + bio10_15_scaled:ScalePET +  
                SnowMonths_cat:ScalePET + scaleAridity:ScalePET + scaleCLYPPT:bio10_15_scaled +  
                scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET + ESA + ScaleElevation +  
                (1 | file/Study_Name), data = bankData, 
              control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))
  
  testData$Predicted <- (predict(mod, testData,  re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$logAbundance, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)
plot(df$predicted ~ df$observed)
abline(0, 1) 


#################################
## BIOMASS
##################################
biomassData <- biomass_model@frame

########
# K-Fold Cross validation
########

k_fold <- 10
splits <- createSplits(biomassData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  
  rows <- as.vector(unlist(splits[k]))
  testData <- biomassData[rows,]
  bankData <- biomassData[-rows,]
  
  mod <- lmer(logBiomass ~ scalePH + scaleCLYPPT + scaleSLTPPT + scaleORCDRC +  
                scaleCECSOL + bio10_7_scaled + bio10_12_scaled + bio10_15_scaled +  
                ScalePET + SnowMonths_cat + scalePH:scaleCLYPPT + scalePH:scaleSLTPPT +  
                scalePH:scaleORCDRC + scalePH:scaleCECSOL + scaleCLYPPT:scaleCECSOL +  
                scaleORCDRC:scaleCECSOL + bio10_7_scaled:bio10_12_scaled +  
                bio10_12_scaled:bio10_15_scaled + bio10_12_scaled:ScalePET +  
                bio10_12_scaled:SnowMonths_cat + bio10_15_scaled:ScalePET +  
                bio10_15_scaled:SnowMonths_cat + ScalePET:SnowMonths_cat +  
                scaleSLTPPT:bio10_12_scaled + scaleCLYPPT:bio10_15_scaled +  
                ESA +
                (1 | file/Study_Name), data = bankData, 
              control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))
  
  testData$Predicted <- (predict(mod, testData, re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$logBiomass, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)
plot(df$predicted ~ df$observed)
abline(0, 1) 


