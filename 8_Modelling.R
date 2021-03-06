########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\USers\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
}


#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
library(lme4)
library(glmmTMB) # github version
library(glmmADMB) # R-forge version
library(car)
library(DHARMa)
source("Functions/FormatData.R")
source("Functions/lme4_ModellingFunctions.R")
source("Functions/ModelSimplification.R")
source("MEE3_1_sm_Appendix_S1/HighstatLib.R")
source("Functions/CorvifVariablePicker.R")

#################################################
# 2. Loading in variables
#################################################

data_in <-"7_Data"

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

if(!dir.exists("8_Data")){
  dir.create("8_Data")
}

data_out <- "8_Data"

if(!dir.exists("Models")){
  dir.create("Models")
}
models <- "Models"
#################################################
# 3. Load in data
#################################################

 sites <- read.csv(file.path(data_in, loadin))
# sites <- read.csv("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\3_Data\\Sites_2017-11-09.csv")
rm(loadin)

#################################################
# 4. Set reference levels
#################################################

sites <- SiteLevels(sites) ## relevels all land use/habitat variables



#################################################
# 4. Species Richness
#################################################

richness <- sites[complete.cases(sites$SpeciesRichness),] #6089 # 7789
richness <- droplevels(richness[richness$ESA != "Unknown",]) # 
richness <- droplevels(richness[-which(richness$SpeciesRichness != round(richness$SpeciesRichness)),]) # 5642 # 7539

# richness <- richness[complete.cases(richness$scalePH),]

# richness <- droplevels(richness[!(is.na(richness$PHIHOX)),])
richness <- droplevels(richness[!(is.na(richness$bio10_15)),]) ## 
richness <- droplevels(richness[!(is.na(richness$OCFinal)),]) ## 
richness <- droplevels(richness[!(is.na(richness$phFinal)),]) ## 
richness <- droplevels(richness[!(is.na(richness$scaleAridity)),]) ## 
richness <- droplevels(richness[!(is.na(richness$SnowMonths_cat)),]) ## 5799 # 7419


table(richness$ESA)
richness_notinclude <- c("Needleleaf deciduous forest", "Tree open",
                         "Sparse vegetation",  "Cropland/Other vegetation mosaic", "Urban",
                         "Bare area (consolidated", "Paddy field", "Wetland/Herbaceous", "Water bodies")

richness <- droplevels(richness[!(richness$ESA %in% richness_notinclude),]) ##    5737 #  7386
summary(richness$phFinal)
richness$scalePH <- as.vector(scale(richness$phFinal))
richness$scaleCLYPPT <- scale(richness$ClayFinal)
richness$scaleSLTPPT <- scale(richness$SiltFinal)
richness$scaleCECSOL <- scale(richness$CECSOL)
richness$scaleORCDRC <- scale(richness$OCFinal)

richness$bio10_1_scaled <- scale(richness$bio10_1)
richness$bio10_4_scaled <- scale(richness$bio10_4)
richness$bio10_7_scaled <- scale(richness$bio10_7)
richness$bio10_12_scaled <- scale(richness$bio10_12)
richness$bio10_15_scaled <- scale(richness$bio10_15)

richness$scaleAridity <- scale(richness$Aridity)
richness$ScalePET <- scale(richness$PETyr)
richness$ScalePETSD <- scale(richness$PET_SD)
richness$scaleElevation <- scale(richness$elevation)

## Save the data
write.csv(richness, file = file.path(data_out, paste("sitesRichness_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

## 
ind <- df_variables(richness)
dat <- richness[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# bio10_7,bio10_15,CECSOL,elevation,Aridity,PETyr,phFinal,ClayFinal,SiltFinal,OCFinal
# same during correction

r1_hurdle <- glmmTMB(SpeciesRichness ~  ESA + scaleElevation + 
                       (scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
                       (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                          ScalePET)^2 + 
                       scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                       scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
                       scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
                       (1|file/Study_Name), data = richness,  family = truncated_poisson,# verbose= TRUE,
                     zi = ~ESA + scaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC +
                      bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                     ScalePET,
                     control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max=2e5)))
save(r1_hurdle, file = file.path(models, "richnessmodel_initialmodel_correction.rds"))



simulationOutput_r1_hurdle <- simulateResiduals(fittedModel = r1_hurdle, n = 250)
plot(simulationOutput_r1_hurdle,quantreg = TRUE)
testZeroInflation(simulationOutput_r1_hurdle, plot = TRUE, alternative = "greater")
testDispersion(simulationOutput_r1_hurdle, alternative = "greater", plot = TRUE)

# Run on cluster
# richness_model <- modelSimplificationAIC_glmmTMB(model = r1_hurdle, itermax = 2e5, evalmax=2e5, dat = richness)
load(file.path(models, "richnessmodel_correction.rds"))
summary(richness_model)

#tt <- getME(r1,"theta")
#ll <- getME(r1,"lower")
#min(tt[ll==0])

# derivs1 <- r1@optinfo$derivs
# sc_grad1 <- with(derivs1,solve(Hessian,gradient))
# max(abs(sc_grad1))
# max(pmin(abs(sc_grad1),abs(derivs1$gradient)))


#################################################
# 5. Biomass
#################################################
biomass <- sites[complete.cases(sites$logBiomass),] # 3689 # 4207 
biomass <- droplevels(biomass[biomass$ESA != "Unknown",]) # 

# biomass <- droplevels(biomass[!(is.na(biomass$PHIHOX)),])
biomass <- droplevels(biomass[!(is.na(biomass$bio10_15)),]) ## 
biomass <- droplevels(biomass[!(is.na(biomass$OCFinal)),]) ## 
biomass <- droplevels(biomass[!(is.na(biomass$phFinal)),]) ## 
biomass <- droplevels(biomass[!(is.na(biomass$SnowMonths_cat)),]) ##  
biomass <- droplevels(biomass[!(is.na(biomass$Aridity)),]) ##  4200


table(biomass$ESA)
biomass_notinclude <- c("Tree open", "Sparse vegetation", "Cropland/Other vegetation mosaic",
                        "Urban", "Paddy field")

biomass <- droplevels(biomass[!(biomass$ESA %in% biomass_notinclude),]) ##   3326 # 4167
summary(biomass$phFinal)
biomass$scalePH <- as.vector(scale(biomass$phFinal))
biomass$scaleCLYPPT <- scale(biomass$ClayFinal)
biomass$scaleSLTPPT <- scale(biomass$SiltFinal)
biomass$scaleCECSOL <- scale(biomass$CECSOL)
biomass$scaleORCDRC <- scale(biomass$OCFinal)

biomass$bio10_1_scaled <- scale(biomass$bio10_1)
biomass$bio10_4_scaled <- scale(biomass$bio10_4)
biomass$bio10_7_scaled <- scale(biomass$bio10_7)
biomass$bio10_12_scaled <- scale(biomass$bio10_12)
biomass$bio10_15_scaled <- scale(biomass$bio10_15)


biomass$scaleAridity <- scale(biomass$Aridity)
biomass$ScalePET <- scale(biomass$PETyr)
biomass$ScalePETSD <- scale(biomass$PET_SD)
biomass$ScaleElevation <- scale(biomass$elevation)


## Save the data
write.csv(biomass, file = file.path(data_out, paste("sitesBiomass_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

## 
ind <- df_variables(biomass)
dat <- biomass[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# bio10_7, bio10_12  bio10_15, CECSOL, 
# elevation, PETyr,phFinal,ClayFinal,SiltFinal,OCFinal

## All fine
# Same variables during correciton phase


## Try with glmmTMB
b1 <- glmmTMB(logBiomass ~  ESA + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
             (bio10_7_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePET + SnowMonths_cat)^2 + 
             scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
             (1|file/Study_Name), data = biomass, 
             ziformula = ~ESA + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL +
               bio10_7_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePET + SnowMonths_cat,
             control = glmmTMBControl(optCtrl = list(iter.max = 2e5,eval.max=2e5))) #, optimizer ="bobyqa"))
## Checking for 'larger hessians'. Small (~0) numbers give non-postivie Hessian matrix warning
ff <- fixef(b1)$zi
round(exp(c(intercept=unname(ff[1]),ff[-1]+ff[1])),3)

# biomass$prediction <- predict(b1)
# plot(biomass$logBiomass ~ biomass$prediction)


simulationOutput_bm <- simulateResiduals(fittedModel = b1, n = 250)
plot(simulationOutput_bm,quantreg = TRUE)
# It seems you are diagnosing a glmmTMB model. There are still a few minor limitations associated with this package. 
# The most important is that glmmTMB doesn't implement an option to create unconditional predictions from the model, 
# which means that predicted values (in res ~ pred) plots include the random effects. With strong random effects, 
# this can sometimes create diagonal patterns from bottom left to top right in the res ~ pred plot. 
# All other tests and plots should work as desired. 
# Please see https://github.com/florianhartig/DHARMa/issues/16 for further details.

# so ignore the res versus pred plot, but look at others
testDispersion(simulationOutput_bm, alternative = "greater", plot = TRUE)
# A bit overdispersed
testZeroInflation(simulationOutput_bm, plot = TRUE, alternative = "greater")
## But zeroinflated


biomass_model <- modelSimplificationAIC_glmmTMB(model = b1, itermax = 2e5, evalmax=2e5, dat = biomass)
  
# biomass_model <- modelSimplificationAIC(model = b1, data = biomass, optimizer = "bobyqa", Iters = 2e5)
# save(biomass_model, file = file.path(models, "biomassmodel_full_revised.rds"))
# load(file.path(models, "biomassmodel_full_revised.rds"))
save(biomass_model, file = file.path(models, "biomassmodel_full_correction.rds"))
# load(file.path(models, "biomassmodel_full_correction.rds"))
# 

simulationOutput_bm <- simulateResiduals(fittedModel = biomass_model, n = 250)
plot(simulationOutput_bm,quantreg = TRUE)
# Pattern looks bad, but is because of GlmmTMB
testDispersion(simulationOutput_bm, alternative = "greater", plot = TRUE)
# a bit overdispersed
testZeroInflation(simulationOutput_bm, plot = TRUE, alternative = "greater")
## And a bit zeroinflated, but not as bad as it was

biomass$predicted <- predict(biomass_model, biomass)
plot(biomass$logBiomass, biomass$predicted)
abline(0, 1)

#################################################
# 6. Abundance
#################################################
hist(sites$logAbundance)
abundance <- sites[complete.cases(sites$logAbundance),] # 7111 # 9089
abundance <- droplevels(abundance[abundance$ESA != "Unknown",]) #

# abundance <- droplevels(abundance[!(is.na(abundance$PHIHOX)),])
abundance <- droplevels(abundance[!(is.na(abundance$bio10_15)),]) ##   
abundance <- droplevels(abundance[!(is.na(abundance$OCFinal)),]) ##  
abundance <- droplevels(abundance[!(is.na(abundance$phFinal)),]) ##  
abundance <- droplevels(abundance[!(is.na(abundance$SnowMonths_cat)),]) ##  
abundance <- droplevels(abundance[!(is.na(abundance$Aridity)),]) ##  6456 # 8724


table(abundance$ESA)
abundance_notinclude <- c("Needleleaf deciduous forest", "Tree open", "Sparse vegetation", "Urban", 
                         "Bare area (consolidated", "Bare area (unconsolidated",  "Paddy field", "Wetland/Herbaceous",
                         "Water bodies")

abundance <- droplevels(abundance[!(abundance$ESA %in% abundance_notinclude),]) #    8677


abundance$scalePH <- as.vector(scale(abundance$phFinal))
abundance$scaleCLYPPT <- scale(abundance$ClayFinal)
abundance$scaleSLTPPT <- scale(abundance$SiltFinal)
abundance$scaleCECSOL <- scale(abundance$CECSOL)
abundance$scaleORCDRC <- scale(abundance$OCFinal)

abundance$bio10_1_scaled <- scale(abundance$bio10_1)
abundance$bio10_4_scaled <- scale(abundance$bio10_4)
abundance$bio10_7_scaled <- scale(abundance$bio10_7)
abundance$bio10_12_scaled <- scale(abundance$bio10_12)
abundance$bio10_15_scaled <- scale(abundance$bio10_15)

abundance$scaleAridity <- scale(abundance$Aridity)
abundance$ScalePET <- scale(abundance$PETyr)
abundance$ScalePETSD <- scale(abundance$PET_SD)
abundance$ScaleElevation  <- scale(abundance$elevation)


## Save the data
write.csv(abundance, file = file.path(data_out, paste("sitesAbundance_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
#abundance <- read.csv(file.path(data_out, "sitesAbundance_2018-09-25.csv"))

ind <- df_variables(abundance)
dat <- abundance[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# bio10_7, bio10_15,CECSOL,elevation,Aridity,PETyr,  phFinal,ClayFinal,SiltFinal, OCFinal   
# same variables as before
a1 <- glmmTMB(logAbundance ~  ESA + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                ScalePET)^2 +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
             scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET + 
             scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity + 
             (1|file/Study_Name), data = abundance, 
             zi = ~ESA + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC +
               bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + ScalePET,
           control = glmmTMBControl(optCtrl = list(iter.max = 2e5,eval.max=2e5)))
# no convergence issues
plot(a1)
simulationOutput_a1 <- simulateResiduals(fittedModel = a1, n = 250)

plot(simulationOutput_a1,quantreg = TRUE)
# simulationOutput_a1b <- simulateResiduals(fittedModel = a1, refit = TRUE)
# ignore these plots
testZeroInflation(simulationOutput_a1, plot = TRUE, alternative = "greater")
## not zeroinflated
testDispersion(simulationOutput_a1, alternative = "greater", plot = TRUE)
# A bit overdispersed


# Run on cluster

abundance_model <- modelSimplificationAIC_glmmTMB(model = a1, itermax = 2e5, evalmax=2e5, dat = abundance)

# load(file.path(models, "abundancemodel_full_correction.rds"))


simulationOutput_a2 <- simulateResiduals(fittedModel = abundance_model, n = 250)
plot(simulationOutput = simulationOutput_a2,quantreg = TRUE)
testZeroInflation(simulationOutput_a2, plot = TRUE, alternative = "greater")
# Nope. none basically

dat <- abundance_model@frame
dat$predicted <- predict(abundance_model, dat)
plot(dat$logAbundance, dat$predicted2)
abline(0, 1)
#########################################################
# Correlation between abundance and biomass
######################################################
#
# As per a reviewers comment
# What is the relationship like between abundance and biomass
# in the raw data.

plot(abundance$logAbundance ~ abundance$logBiomass)

abundancevBiomass <- lmer(logAbundance ~  logBiomass + 
             (1 + logBiomass|file/Study_Name), data = abundance,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))


summary(abundancevBiomass)
r.squaredGLMM(abundancevBiomass)

coef(abundancevBiomass)$Study_Name
