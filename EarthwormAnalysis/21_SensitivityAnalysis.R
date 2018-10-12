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


figures <- "Figures"
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

## FG Richness
data_in <- "15_Data"

files <- list.files(file.path(data_in))
files <- files[grep("sites\\+FGRichness_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
fgrichness <- read.csv(file.path(data_in, loadin))



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
load(file.path(models, "fgrichnessmodel.rds"))



k_fold <- 10



#################################################
# 5. Species Richness
################################################
richnessData <- richness_model@frame


########
# K-Fold Cross validation
########

splits <- createSplits(richnessData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  
  rows <- as.vector(unlist(splits[k]))
  testData <- richnessData[rows,]
  bankData <- richnessData[-rows,]
  
  mod <-  glmer(formula = richness_model@call$formula, data = bankData, family = "poisson",
                control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
  testData$Predicted <- (predict(mod, testData,  re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$SpeciesRichness, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)
plot(df$predicted ~ df$observed)
abline(0, 1) 

richness <- df

#################################################
# 5. Abundance
################################################

abundanceData <- abundance_model@frame


########
# K-Fold Cross validation
########

splits <- createSplits(abundanceData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  
  rows <- as.vector(unlist(splits[k]))
  testData <- abundanceData[rows,]
  bankData <- abundanceData[-rows,]
  
  mod <-  lmer(formula = abundance_model@call$formula, data = bankData, 
               control = lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
  testData$Predicted <- (predict(mod, testData,  re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$logAbundance, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)

# jpeg(file = file.path(figures, "Abundance_Crossvalidation.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
plot(df$predicted ~ df$observed, ylab = "Predicted log-Abundance", xlab = "Observed log-Abundance", pch = 19, cex = 0.5)
abline(0, 1) 
# dev.off()

abundance <- df

#################################
## BIOMASS
##################################
biomassData <- biomass_model@frame

########
# K-Fold Cross validation
########

splits <- createSplits(biomassData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  
  rows <- as.vector(unlist(splits[k]))
  testData <- biomassData[rows,]
  bankData <- biomassData[-rows,]
  
  mod <-  lmer(formula = biomass_model@call$formula, data = bankData, 
               control = lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
           
  testData$Predicted <- (predict(mod, testData, re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$logBiomass, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)
plot(df$predicted ~ df$observed)
abline(0, 1) 

biomass <- df
#################################
## FG Richness
##################################
fgRichnessData <- fgrichness_model@frame

########
# K-Fold Cross validation
########
splits <- createSplits(fgRichnessData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  
  rows <- as.vector(unlist(splits[k]))
  testData <- fgRichnessData[rows,]
  bankData <- fgRichnessData[-rows,]
  
  mod <-  glmer(formula = fgrichness_model@call$formula, data = bankData, family = "poisson",
                control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
  testData$Predicted <- (predict(mod, testData, re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$FGRichness, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)
plot(df$predicted ~ df$observed)
abline(0, 1) 

fgrichness <- df

####################################
## PLOT
####################################

jpeg(file = file.path(figures, "AllModels_Crossvalidation.jpg"), quality = 100, res = 200, width = 2000, height = 1500)

par(mar = c(2.5, 2.5, 1, 1))
par(mfrow = c(2, 2))
plot(exp(richness$predicted) ~ jitter(richness$observed), ylab = "", xlab = "", pch = 19, cex = 0.5)
abline(0, 1) 
text(x = -0.5, y = 12, labels = "Species Richness", pos = 4)
plot(abundance$predicted ~ abundance$observed, ylab = "", xlab = "", pch = 19, cex = 0.5, ylim = c(0, 8))
abline(0, 1) 
text(x = -0.2, y = 7.4, labels = "(log)Abundance", pos = 4)
plot(biomass$predicted ~ biomass$observed, ylab = "", xlab = "", pch = 19, cex = 0.5, ylim = c(0, 8))
abline(0, 1) 
text(x = -0.2, y = 7, labels = "(log)Biomass", pos = 4)

plot(exp(fgrichness$predicted) ~ jitter(fgrichness$observed), ylab = "", xlab = "", pch = 19, cex = 0.5)
abline(0, 1) 
text(x = -0.3, y = 4, labels = "Functional Richness", pos = 4)

dev.off()
