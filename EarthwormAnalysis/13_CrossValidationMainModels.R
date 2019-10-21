if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\Users\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis")
}

#################################################
# 1. Libraries
#################################################
library(lme4)
library(MuMIn)

source(file.path("Functions", "CrossValidationAndMSE.R"))

figures <- "Figures"

#################################################
# 3. Create directories
#################################################

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

if(!dir.exists("13_Data")){
  dir.create("13_Data")
}
data_out <- "13_Data"

#################################################
# 4. Load in models
#################################################

models <- "Models"


load(file.path(models, "richnessmodel_revised.rds"))
load(file.path(models, "biomassmodel_full_revised.rds"))
load(file.path(models, "abundancemodel_full_revised.rds"))


k_fold <- 10


abundanceData <- abundance_model@frame
biomassData <- biomass_model@frame
richnessData <- richness_model@frame

#################################################
# NOT ACCOUNTING STUDIES
################################################

#################################################
# 5. Abundance
################################################


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
optimizer = "bobyqa"
Iters = 2e5
data <- read.csv(file.path("8_Data", "sitesRichness_2019-06-20.csv", sep = ""))
fam = "poisson"
r.squaredGLMM(richness_model)




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

richness$predicted <- exp(richness$predicted) 

calculateMSE(richness)
calculateMSEofQuantiles(richness)

write.csv(richness, file = file.path(data_out, "RichnessCrossValidation.csv"), row.names = FALSE)



####################################
## PLOT
####################################
richness <- read.csv(file.path(data_out, "RichnessCrossValidation.csv"))
abundance <- read.csv(file = file.path(data_out, "AbundanceCrossValidation.csv"))
biomass <- read.csv(file = file.path(data_out, "BiomassCrossValidation.csv"))



jpeg(file = file.path(figures, "AllModels_Crossvalidation.jpg"), quality = 100, res = 200, width = 3000, height = 1000)
par(oma = c(2, 2, 1, 0))

par(mar = c(3.5, 3.5, 2, 1))
par(mfrow = c(1, 3))

plot(richness$predicted ~ jitter(richness$observed), ylab = "", xlab = "", pch = 19, cex = 0.5)
abline(0, 1) 
# text(x = -0.5, y = 12, labels = "Species Richness", pos = 4)
mtext("Species Richness", side = 3, line = 0.5, at = 0, adj = 0.1)
mtext("A", side = 3, line = 0.5, at = 0, adj = 3, font = 2)

plot(log(abundance$predicted + 1) ~ log(abundance$observed + 1), ylab = "", xlab = "", pch = 19, cex = 0.5)
abline(0, 1) 
# text(x = -0.2, y = 6, labels = "(log)Abundance", pos = 4)
mtext("(ln-) Abundance", side = 3, line = 0.5, at = 0, adj = 0.1)
mtext("B", side = 3, line = 0.5, at = 0, adj = 3, font = 2)


plot(log(biomass$predicted+1) ~ log(biomass$observed + 1), ylab = "", xlab = "", pch = 19, cex = 0.5)
abline(0, 1) 
#text(x = -0.2, y = 5, labels = "(log)Biomass", pos = 4)
mtext("(ln-) Biomass", side = 3, line = 0.5, at = 0, adj = 0.1)
mtext("C", side = 3, line = 0.5, at = 0, adj = 3, font = 2)

mtext('Predicted values', side = 2, outer = TRUE, line = -0, las = 0, cex = 1.5)
mtext('Observed values', side = 1, outer = TRUE, line = -0, las = 0, cex = 1.5)

dev.off()

#################################################
# ACCOUNTING FOR STUDIES
################################################


#################################################
# abundance
################################################

splits <- createSplitsStudies(abundanceData, kfold = 10)


predictedData <- list()
for(k in 1:k_fold){
  
  studies <- as.vector(unlist(splits[k]))
  testData <- abundanceData[abundanceData$Study_Name %in% studies,]
  bankData <- abundanceData[!(abundanceData$Study_Name %in% studies),]
  
  mod <-  lmer(formula = abundance_model@call$formula, data = bankData, 
               control = lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
  
  testData$Predicted <- (predict(mod, testData, re.form = NULL, allow.new.levels = TRUE))
  
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

write.csv(abundance, file = file.path(data_out, "AbundanceStudyCrossValidation.csv"), row.names = FALSE)

#################################################
# biomass
################################################

splits <- createSplitsStudies(biomassData, kfold = 10)

biomassData <- biomassData[biomassData$ESA != "Bare area (unconsolidated",] ## Actually not included in the model anyway

predictedData <- list()
for(k in 1:k_fold){
  
  studies <- as.vector(unlist(splits[k]))
  testData <- biomassData[biomassData$Study_Name %in% studies,]
  bankData <- biomassData[!(biomassData$Study_Name %in% studies),]
  
  mod <-  lmer(formula = biomass_model@call$formula, data = bankData, 
               control = lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
  
  testData$Predicted <- (predict(mod, testData, re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$logBiomass, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)

# jpeg(file = file.path(figures, "Abundance_Crossvalidation.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
plot(df$predicted ~ df$observed, ylab = "Predicted log-Biomass", xlab = "Observed log-Biomass", pch = 19, cex = 0.5)
abline(0, 1) 
# dev.off()

biomass <- df

biomass$predicted <- exp(biomass$predicted) - 1
biomass$observed <- exp(biomass$observed) - 1

calculateMSE(biomass)
calculateMSEofQuantiles(biomass)

write.csv(biomass, file = file.path(data_out, "BiomassStudyCrossValidation.csv"), row.names = FALSE)

#############################################
# richness
############################################


splits <- createSplitsStudies(richnessData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  print(k)
  
  studies <- as.vector(unlist(splits[k]))
  testData <- richnessData[richnessData$Study_Name %in% studies,]
  bankData <- richnessData[!(richnessData$Study_Name %in% studies),]
  
  mod <-  glmer(formula = richness_model@call$formula, data = bankData, family = "poisson",
                control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
  testData$Predicted <- (predict(mod, testData,  re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$SpeciesRichness, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)
plot(df$predicted ~ df$observed)
abline(0, 1) 

richness <- df

richness$predicted <- exp(richness$predicted) 

calculateMSE(richness)
calculateMSEofQuantiles(richness)

write.csv(richness, file = file.path(data_out, "RichnessStudyCrossValidation.csv"), row.names = FALSE)



richness <- read.csv(file.path(data_out, "RichnessStudyCrossValidation.csv"))
abundance <- read.csv(file = file.path(data_out, "AbundanceStudyCrossValidation.csv"))
biomass <- read.csv(file = file.path(data_out, "BiomassStudyCrossValidation.csv"))



jpeg(file = file.path(figures, "AllModels_StudyCrossvalidation.jpg"), quality = 100, res = 200, width = 3000, height = 1000)
par(oma = c(2, 2, 1, 0))

par(mar = c(3.5, 3.5, 2, 1))
par(mfrow = c(1, 3))

plot(richness$predicted ~ jitter(richness$observed), ylab = "", xlab = "", pch = 19, cex = 0.5)
abline(0, 1) 
# text(x = -0.5, y = 12, labels = "Species Richness", pos = 4)
mtext("Species Richness", side = 3, line = 0.5, at = 0, adj = 0.1)
mtext("A", side = 3, line = 0.5, at = 0, adj = 3, font = 2)

plot(log(abundance$predicted + 1) ~ log(abundance$observed + 1), ylab = "", xlab = "", pch = 19, cex = 0.5)
abline(0, 1) 
# text(x = -0.2, y = 6, labels = "(log)Abundance", pos = 4)
mtext("(ln-) Abundance", side = 3, line = 0.5, at = 0, adj = 0.1)
mtext("B", side = 3, line = 0.5, at = 0, adj = 3, font = 2)

plot(log(biomass$predicted+1) ~ log(biomass$observed + 1), ylab = "", xlab = "", pch = 19, cex = 0.5)
abline(0, 1) 
#text(x = -0.2, y = 5, labels = "(log)Biomass", pos = 4)
mtext("(ln-) Biomass", side = 3, line = 0.5, at = 0, adj = 0.1)
mtext("C", side = 3, line = 0.5, at = 0, adj = 3, font = 2)

mtext('Predicted values', side = 2, outer = TRUE, line = -0, las = 0, cex = 1.5)
mtext('Observed values', side = 1, outer = TRUE, line = -0, las = 0, cex = 1.5)

dev.off()




#################################################
# ACCOUNTING FOR STUDIES WITHIN BIOMES
################################################


data_in <- "8_Data"

#################################################
# abundance
################################################

files <- list.files(file.path(data_in))
files <- files[grep('sitesAbundance_', files)]

file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

abundance <- read.csv(file.path(data_in, loadin))
###

# First check that all studies are only in one biome
dfbiome <- data.frame(tapply(abundance$biome, abundance$Study_Name, FUN = function(x) length(unique(x))))
dfbiome <- dfbiome[dfbiome$tapply.abundance.biome..abundance.Study_Name..FUN...function.x..length.unique.x... > 1,]
# 23 studies cross biomes

abundance_tmp <- droplevels(abundance[abundance$Study_Name %in% rownames(dfbiome),])
table(abundance_tmp$Study_Name, abundance_tmp$biome)

### What we can do, is use the whole study regardless of whether some sites are not in the biome

dfbiome <- data.frame(tapply(abundance$Study_Name, abundance$biome, FUN = function(x) length(unique(x))))
abundance$region <- as.character(abundance$country)

europe <- c("Austria", "Belgium", "Croatia", "Czech Republic", "Denmark", "Estonia", "Finland",
            "France", "Germany", "Greece", "Hungary", "Ireland", "Italy", "Lithuania", 'Netherlands',
            "Norway", "Poland", "Portugal", "Romania", "Russia", "Slovakia", "Slovenia", "Spain",
            "Sweden", "Switzerland", "United Kingdom")
abundance$region[abundance$region %in% europe] <- "Europe"
africa <- c("Benin", "Burkina Faso", "Cameroon", "Cote d'Ivoire", "Ghana", "Kenya", "Madagascar", 
            "Malawi", "Nigeria")
abundance$region[abundance$region %in% africa] <- "Africa"
asia <- c("China", "India", "Iran (Islamic Republic of)", "Japan", "Libyan Arab Jamahiriya", "Malaysia", 
          "Philippines")
abundance$region[abundance$region %in% asia] <- "Asia"

northA <- c("Canada", "United States")
abundance$region[abundance$region %in% northA] <- "North America"

southA <- c("Argentina", "Bolivia", "Brazil", "Colombia", "Ecuador", "Honduras", "Mexico", "Nicaragua",
            "Peru", "Puerto Rico", "Uruguay", "Venezuela")
abundance$region[abundance$region %in% southA] <- "South America"


# Australiasia <- c("Australia")

europe <- droplevels(abundance[which(abundance$region == "Europe"),])
splits <- createSplitsStudies(europe, kfold = 10)


predictedData <- list()
for(k in 1:k_fold){
  
  studies <- as.vector(unlist(splits[k]))
  testData <- europe[europe$Study_Name %in% studies,]
  bankData <- europe[!(europe$Study_Name %in% studies),]
  
  mod <-  lmer(formula = abundance_model@call$formula, data = bankData, 
               control = lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
  
  testData$Predicted <- (predict(mod, testData, re.form = NULL, allow.new.levels = TRUE))
  
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

write.csv(abundance, file = file.path(data_out, "AbundanceStudyCrossValidation.csv"), row.names = FALSE)

#################################################
# biomass
################################################

splits <- createSplitsStudies(biomassData, kfold = 10)

biomassData <- biomassData[biomassData$ESA != "Bare area (unconsolidated",] ## Actually not included in the model anyway

predictedData <- list()
for(k in 1:k_fold){
  
  studies <- as.vector(unlist(splits[k]))
  testData <- biomassData[biomassData$Study_Name %in% studies,]
  bankData <- biomassData[!(biomassData$Study_Name %in% studies),]
  
  mod <-  lmer(formula = biomass_model@call$formula, data = bankData, 
               control = lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
  
  testData$Predicted <- (predict(mod, testData, re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$logBiomass, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)

# jpeg(file = file.path(figures, "Abundance_Crossvalidation.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
plot(df$predicted ~ df$observed, ylab = "Predicted log-Biomass", xlab = "Observed log-Biomass", pch = 19, cex = 0.5)
abline(0, 1) 
# dev.off()

biomass <- df

biomass$predicted <- exp(biomass$predicted) - 1
biomass$observed <- exp(biomass$observed) - 1

calculateMSE(biomass)
calculateMSEofQuantiles(biomass)

write.csv(biomass, file = file.path(data_out, "BiomassStudyCrossValidation.csv"), row.names = FALSE)

#############################################
# richness
############################################


splits <- createSplitsStudies(richnessData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  print(k)
  
  studies <- as.vector(unlist(splits[k]))
  testData <- richnessData[richnessData$Study_Name %in% studies,]
  bankData <- richnessData[!(richnessData$Study_Name %in% studies),]
  
  mod <-  glmer(formula = richness_model@call$formula, data = bankData, family = "poisson",
                control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
  
  testData$Predicted <- (predict(mod, testData,  re.form = NULL, allow.new.levels = TRUE))
  
  predictedData[[k]] <- data.frame(observed = testData$SpeciesRichness, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)
plot(df$predicted ~ df$observed)
abline(0, 1) 

richness <- df

richness$predicted <- exp(richness$predicted) 

calculateMSE(richness)
calculateMSEofQuantiles(richness)

write.csv(richness, file = file.path(data_out, "RichnessStudyCrossValidation.csv"), row.names = FALSE)



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
