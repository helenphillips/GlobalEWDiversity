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
source("Functions/CorvifVariablePicker.R")

#################################################
# 2. Loading in variables
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

#################################################
# 2.5 Create folders
#################################################

if(!dir.exists("4_Data")){
  dir.create("4_Data")
}

data_out <- "4_Data"

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

richness <- sites[complete.cases(sites$SpeciesRichness),] #6089
richness <- droplevels(richness[richness$ESA != "Unknown",]) # 5660
richness <- droplevels(richness[-which(richness$SpeciesRichness != round(richness$SpeciesRichness)),]) # 5642

# richness <- richness[complete.cases(richness$scalePH),]

# richness <- droplevels(richness[!(is.na(richness$PHIHOX)),])
richness <- droplevels(richness[!(is.na(richness$bio10_15)),]) ## 5629
richness <- droplevels(richness[!(is.na(richness$OCFinal)),]) ## 5622
richness <- droplevels(richness[!(is.na(richness$phFinal)),]) ## 5618
richness <- droplevels(richness[!(is.na(richness$scaleAridity)),]) ## 5509
richness <- droplevels(richness[!(is.na(richness$SnowMonths_cat)),]) ## 5466


table(richness$ESA)
richness_notinclude <- c("Needleleaf deciduous forest", "Tree open",
                         "Sparse vegetation",  "Cropland/Other vegetation mosaic",
                         "Bare area (consolidated", "Paddy field", "Wetland/Herbaceous", "Water bodies")

richness <- droplevels(richness[!(richness$ESA %in% richness_notinclude),]) ##   5414
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

r1 <- glmer(SpeciesRichness ~  ESA + scaleElevation + (scalePH  + 
             scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                ScalePET)^2 + 
              scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
              scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
              scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
             (1|file/Study_Name), data = richness, family = poisson,
            control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))


summary(r1)
save(r1, file = file.path(models, "richnessmodel_initialmodel.rds"))
load(file.path(models, "richnessmodel_initialmodel.rds"))


simulationOutput_r1 <- simulateResiduals(fittedModel = r1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput_r1,quantreg = TRUE)
# simulationOutput_r1b <- simulateResiduals(fittedModel = r1, refit = TRUE)
testDispersion(simulationOutput_r1, alternative = "greater", plot = TRUE)
testZeroInflation(simulationOutput_r1, plot = TRUE, alternative = "greater")
## No zero inflation or overdispersion

#tt <- getME(r1,"theta")
#ll <- getME(r1,"lower")
#min(tt[ll==0])

# derivs1 <- r1@optinfo$derivs
# sc_grad1 <- with(derivs1,solve(Hessian,gradient))
# max(abs(sc_grad1))
# max(pmin(abs(sc_grad1),abs(derivs1$gradient)))

richness_model <- modelSimplificationAIC(model = r1, data = richness, optimizer = "bobyqa", Iters = 2e5)
save(richness_model, file = file.path(models, "richnessmodel_full.rds"))
# load(file.path(models, "richnessmodel_full.rds"))

##  From DHARMa
simulationOutput <- simulateResiduals(fittedModel = richness_model, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput,quantreg = TRUE)

#################################################
# 5. Biomass
#################################################
biomass <- sites[complete.cases(sites$logBiomass),] # 3689
biomass <- droplevels(biomass[biomass$ESA != "Unknown",]) # 3368

# biomass <- droplevels(biomass[!(is.na(biomass$PHIHOX)),])
biomass <- droplevels(biomass[!(is.na(biomass$bio10_15)),]) ## 3365
biomass <- droplevels(biomass[!(is.na(biomass$OCFinal)),]) ## 3364
biomass <- droplevels(biomass[!(is.na(biomass$phFinal)),]) ## 3364
biomass <- droplevels(biomass[!(is.na(biomass$SnowMonths_cat)),]) ##  3361
biomass <- droplevels(biomass[!(is.na(biomass$Aridity)),]) ##  3357


table(biomass$ESA)
biomass_notinclude <- c("Tree open", "Sparse vegetation", "Cropland/Other vegetation mosaic",
                        "Urban", "Paddy field")

biomass <- droplevels(biomass[!(biomass$ESA %in% biomass_notinclude),]) ##   3324
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

b1 <- lmer(logBiomass ~  ESA + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                (bio10_7_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePET + SnowMonths_cat)^2 + 
            scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                (1|file/Study_Name), data = biomass,
              control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

simulationOutput <- simulateResiduals(fittedModel = b1, n = 250)
plot(simulationOutput,quantreg = TRUE)
# pretty good

biomass_model <- modelSimplificationAIC(model = b1, data = biomass, optimizer = "bobyqa", Iters = 2e5)
save(biomass_model, file = file.path(models, "biomassmodel_full.rds"))
# load(file.path(models, "biomassmodel_full.rds"))

simulationOutput_bm <- simulateResiduals(fittedModel = biomass_model, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput_bm,quantreg = TRUE)
testDispersion(simulationOutput_bm, alternative = "greater", plot = FALSE)
# Not overdispersed
testZeroInflation(simulationOutput_bm, plot = TRUE, alternative = "greater")
## But zeroinflated

dat <- biomass_model@frame
dat$predicted <- predict(biomass_model, dat)
plot(dat$logBiomass, dat$predicted)
abline(0, 1)

#################################################
# 6. Abundance
#################################################
hist(sites$logAbundance)
abundance <- sites[complete.cases(sites$logAbundance),] # 7211
abundance <- droplevels(abundance[abundance$ESA != "Unknown",]) #6759

# abundance <- droplevels(abundance[!(is.na(abundance$PHIHOX)),])
abundance <- droplevels(abundance[!(is.na(abundance$bio10_15)),]) ##   
abundance <- droplevels(abundance[!(is.na(abundance$OCFinal)),]) ##  
abundance <- droplevels(abundance[!(is.na(abundance$phFinal)),]) ##  6731
abundance <- droplevels(abundance[!(is.na(abundance$SnowMonths_cat)),]) ##  6657
abundance <- droplevels(abundance[!(is.na(abundance$Aridity)),]) ##  6576


table(abundance$ESA)
abundance_notinclude <- c("Needleleaf deciduous forest", "Tree open", "Sparse vegetation", 
                         "Bare area (consolidated", "Bare area (unconsolidated",  "Paddy field", "Wetland/Herbaceous",
                         "Water bodies")

abundance <- droplevels(abundance[!(abundance$ESA %in% abundance_notinclude),]) #  6709


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

ind <- df_variables(abundance)
dat <- abundance[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# bio10_7, bio10_15,CECSOL,elevation,Aridity,PETyr,  phFinal,ClayFinal,SiltFinal, OCFinal   

a1 <- lmer(logAbundance ~  ESA + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                ScalePET)^2 +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
             scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET + 
             scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity + 
             (1|file/Study_Name), data = abundance,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

plot(a1)
simulationOutput_a1 <- simulateResiduals(fittedModel = a1, n = 250)
plot(simulationOutput = simulationOutput_a1,quantreg = TRUE)
# simulationOutput_a1b <- simulateResiduals(fittedModel = a1, refit = TRUE)

testZeroInflation(simulationOutput_a1, plot = TRUE, alternative = "greater")
## But zeroinflated


abundance_model <- modelSimplificationAIC(model = a1, data = abundance, optimizer = "bobyqa", Iters = 2e5)
save(abundance_model, file = file.path(models, "abundancemodel_full.rds"))

simulationOutput_a2 <- simulateResiduals(fittedModel = abundance_model, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput_a2,quantreg = TRUE)
testZeroInflation(simulationOutput_a2, plot = FALSE, alternative = "more")
## But zeroinflated

dat <- abundance_model@frame
dat$predicted <- predict(abundance_model, dat)
plot(dat$logAbundance, dat$predicted)
abline(0, 1)