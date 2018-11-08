if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}
#################################################
# 1. Libraries
#################################################
library(lme4)
source(file.path("Functions", "CrossValidationAndMSE.R"))

figures <- "Figures"

#################################################
# 3. Create directories
#################################################

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

if(!dir.exists("21_Data")){
  dir.create("21_Data")
}
data_out <- "21_Data"

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
# 5. Abundance
################################################

abundanceData <- abundance_model@frame
r.squaredGLMM(abundance_model)


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

abundance$predicted <- exp(abundance$predicted) - 1
abundance$observed <- exp(abundance$observed) - 1

calculateMSE(abundance)
calculateMSEofQuantiles(abundance)

write.csv(abundance, file = file.path(data_out, "AbundanceCrossValidation.csv"), row.names = FALSE)

#################################
## BIOMASS
##################################
r.squaredGLMM(biomass_model)

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


biomass$observed <- exp(biomass$observed) - 1
biomass$predicted <- exp(biomass$predicted) - 1

calculateMSE(biomass)
calculateMSEofQuantiles(biomass)

write.csv(biomass, file = file.path(data_out, "BiomassCrossValidation.csv"), row.names = FALSE)

####
##
## WARNING: THESE TAKE A LONG TIME
##
###
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

richness$observed <- exp(richness$observed)
richness$predicted <- exp(richness$predicted) 

calculateMSE(richness)
calculateMSEofQuantiles(richness)

write.csv(richness, file = file.path(data_out, "RichnessCrossValidation.csv"), row.names = FALSE)

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

fgrichness$observed <- exp(fgrichness$observed)
fgrichness$predicted <- exp(fgrichness$predicted) 

calculateMSE(fgrichness)
calculateMSEofQuantiles(fgrichness)

write.csv(fgrichness, file = file.path(data_out, "FGRichnessCrossValidation.csv"), row.names = FALSE)


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


##################################################
# CHECKING A CHANGE IN CLIMATE WITHIN EACH STUDY
##################################################

allStudies <- c(as.vector(abundance_model@frame$Study_Name), as.vector(biomass_model@frame$Study_Name), 
                as.vector(richness_model@frame$Study_Name), as.vector(fgrichness_model@frame$Study_Name))
allStudies <- unique(allStudies)


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

usedSites <- sites[sites$Study_Name %in% allStudies,]

library(plyr)
library(dplyr)

var.df <- usedSites %>% # Start by defining the original dataframe, AND THEN...
  group_by(Study_Name) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    bio10_1 = var(bio10_1),
    bio10_4 = var(bio10_4),
    bio10_7 = var(bio10_7),
    bio10_12 = var(bio10_12),
    bio10_15 = var(bio10_15)
  )

vardf <- as.data.frame(var.df)

someClimate <- vardf[apply(vardf [c('bio10_1','bio10_4','bio10_7', 'bio10_12', 'bio10_15')],1,function(x) any(x == 0)),]
noClimate <- vardf[apply(vardf [c('bio10_1','bio10_4','bio10_7', 'bio10_12', 'bio10_15')],1,function(x) all(x == 0)),]
