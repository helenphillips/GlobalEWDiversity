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
# 5. Investigating Climate Variables
#################################################
## No longer doing random forest models

## Bio 1, 4, 7, 12, 15

## Checking for collinearity
# plot(sites$bio10_14 ~ sites$bio10_13)
# r1 <- lm(sites$bio10_14 ~ sites$bio10_13)
# summary(r1) # low r1, with no random effect structure

##https://stats.stackexchange.com/questions/82984/how-to-test-and-avoid-multicollinearity-in-mixed-linear-model/142032
## That page said a cut of 4. We are not higher than 0.9
x <- data.frame(sites$bio10_1,sites$bio10_4,sites$bio10_7,sites$bio10_12,sites$bio10_15, 
                sites$SnowMonths, sites$Aridity, sites$PETyr, sites$PET_SD)
correl_dummy_df <- round(cor(x, use = "pair"), 2) 

## 7 is highly collinear with bio4 and quite collinear with bio1

## VIFs
source("MEE3_1_sm_Appendix_S1/HighstatLib.R")
corvif(x) ## There might be an issue here
y <- data.frame(sites$bio10_1,sites$bio10_4,sites$bio10_12,sites$bio10_15,
                sites$SnowMonths, sites$Aridity, sites$PETyr, sites$PET_SD)
corvif(y)
corvif(x)
## Problems with bio4 and bio7

# Maybe problem with PET
z <- data.frame(sites$bio10_1,sites$bio10_4,sites$bio10_12,sites$bio10_15,
                sites$SnowMonths, sites$Aridity, sites$PET_SD)
corvif(z)

aa <-  data.frame(sites$bio10_1,sites$bio10_4,sites$bio10_12,sites$bio10_15,
                  sites$Aridity, sites$PETyr, sites$PET_SD)

corvif(aa)

bb <- data.frame(sites$bio10_1,sites$bio10_4,sites$bio10_12,sites$bio10_15,
                sites$SnowMonths, sites$Aridity)
corvif(bb)

cc <- data.frame(sites$bio10_1,sites$bio10_4,sites$bio10_12,sites$bio10_15,
                 sites$SnowMonths)
corvif(cc)

dd <- data.frame(sites$bio10_1,sites$bio10_4,sites$bio10_12,sites$bio10_15,
                 sites$Aridity)
corvif(dd)


#################################################
# 4. Species Richness
#################################################

richness <- sites[complete.cases(sites$SpeciesRichness),] #6093
richness <- droplevels(richness[richness$ESA != "Unknown",]) # 5664
richness <- droplevels(richness[-which(richness$SpeciesRichness != round(richness$SpeciesRichness)),]) # 5646

# richness <- richness[complete.cases(richness$scalePH),]

# richness <- droplevels(richness[!(is.na(richness$PHIHOX)),])
richness <- droplevels(richness[!(is.na(richness$bio10_1)),]) ## 5633
richness <- droplevels(richness[!(is.na(richness$OCFinal)),]) ## 5622
richness <- droplevels(richness[!(is.na(richness$phFinal)),]) ## 5622


table(richness$ESA)
richness_notinclude <- c("Needleleaf deciduous forest", "Tree open",
                         "Sparse vegetation",  "Cropland/Other vegetation mosaic",
                         "Bare area (consolidated", "Paddy field", "Wetland/Herbaceous", "Water bodies")

richness <- droplevels(richness[!(richness$ESA %in% richness_notinclude),]) ##  5565
summary(richness$phFinal)
richness$scalePH <- as.vector(scale(richness$phFinal))
richness$scaleCLYPPT <- scale(richness$ClayFinal)
richness$scaleSLTPPT <- scale(richness$SiltFinal)
richness$scaleCECSOL <- scale(richness$CECSOL)
richness$scaleORCDRC <- scale(richness$OCFinal)


# Three most important climate variables
richness$bio10_1_scaled <- scale(richness$bio10_1)
richness$bio10_4_scaled <- scale(richness$bio10_4)
# richness$bio10_7_scaled <- scale(richness$bio10_7)
richness$bio10_12_scaled <- scale(richness$bio10_12)
richness$bio10_15_scaled <- scale(richness$bio10_15)

richness$scaleAridity <- scale(richness$Aridity)
richness$ScalePET <- scale(richness$PETyr)
richness$ScalePETSD <- scale(richness$PET_SD)

## Save the data
write.csv(richness, file = file.path(data_out, paste("sitesRichness_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

r1 <- glmer(SpeciesRichness ~  ESA + (scalePH  + 
             scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_1_scaled + bio10_4_scaled + 
                bio10_12_scaled + bio10_15_scaled +
                SnowMonths)^2 + 
             #  SNDPPT # Not included, as the other two dictate the third
              
             #  (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = richness, family = poisson,
            control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))


summary(r1)

simulationOutput_r1 <- simulateResiduals(fittedModel = r1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput_r1,quantreg = TRUE)
# simulationOutput_r1b <- simulateResiduals(fittedModel = r1, refit = TRUE)
testOverdispersion(simulationOutput_r1, alternative = "overdispersion", plot = TRUE)
testZeroInflation(simulationOutput_r1, plot = TRUE, alternative = "more")
## No zero inflation or overdispersion

#tt <- getME(r1,"theta")
#ll <- getME(r1,"lower")
#min(tt[ll==0])

# derivs1 <- r1@optinfo$derivs
# sc_grad1 <- with(derivs1,solve(Hessian,gradient))
# max(abs(sc_grad1))
# max(pmin(abs(sc_grad1),abs(derivs1$gradient)))

richness_model <- modelSimplification(model = r1, data = richness, alpha = 0.05, optimizer = "bobyqa", Iters = 2e5)
save(richness_model, file = file.path(models, "richnessmodel_full.rds"))
# load(file.path(models, "richnessmodel_full.rds"))



##  From DHARMa
simulationOutput <- simulateResiduals(fittedModel = richness_model, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput,quantreg = TRUE)

#################################################
# 5. Biomass
#################################################
biomass <- sites[complete.cases(sites$logBiomass),] # 3698
biomass <- droplevels(biomass[biomass$ESA != "Unknown",]) # 3377

# biomass <- droplevels(biomass[!(is.na(biomass$PHIHOX)),])
biomass <- droplevels(biomass[!(is.na(biomass$bio10_1)),]) ## 3374
biomass <- droplevels(biomass[!(is.na(biomass$OCFinal)),]) ## 3373
biomass <- droplevels(biomass[!(is.na(biomass$phFinal)),]) ## 3373


table(biomass$ESA)
biomass_notinclude <- c("Tree open", "Sparse vegetation", "Cropland/Other vegetation mosaic",
                        "Urban", "Paddy field")

biomass <- droplevels(biomass[!(biomass$ESA %in% biomass_notinclude),]) ##   3339
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


## Save the data
write.csv(biomass, file = file.path(data_out, paste("sitesBiomass_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)



b1 <- lmer(logBiomass ~  ESA + (scalePH  + 
            scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_1_scaled + bio10_4_scaled + 
                bio10_12_scaled + bio10_15_scaled + 
                SnowMonths)^2 +
             #  SNDPPT # Not included, as the other two dictate the third
             
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = biomass,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))
plot(b1)

simulationOutput <- simulateResiduals(fittedModel = b1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput,quantreg = TRUE)

biomass_model <- modelSimplification(model = b1, data = biomass, alpha = 0.05, optimizer = "bobyqa", Iters = 2e5)
save(biomass_model, file = file.path(models, "biomassmodel_full.rds"))
# load(file.path(models, "biomassmodel_full.rds"))

simulationOutput_bm <- simulateResiduals(fittedModel = biomass_model, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput_bm,quantreg = TRUE)
testOverdispersion(simulationOutput_bm, alternative = "overdispersion", plot = TRUE)
# Not overdispersed
testZeroInflation(simulationOutput_bm, plot = TRUE, alternative = "more")
## But zeroinflated


#################################################
# 6. Abundance
#################################################
hist(sites$Site_Abundance)
hist(sites$logAbundance)
abundance <- sites[complete.cases(sites$Site_Abundance),]
abundance <- droplevels(abundance[abundance$ESA != "Unknown",]) #6811

# abundance <- droplevels(abundance[!(is.na(abundance$PHIHOX)),])
abundance <- droplevels(abundance[!(is.na(abundance$bio10_1)),]) ##   6796
abundance <- droplevels(abundance[!(is.na(abundance$OCFinal)),]) ##  6783
abundance <- droplevels(abundance[!(is.na(abundance$phFinal)),]) ##  6783


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

## Save the data
write.csv(abundance, file = file.path(data_out, paste("sitesAbundance_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)


a1 <- lmer(logAbundance ~  ESA + (scalePH  + 
                                    scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_1_scaled + bio10_4_scaled + 
                bio10_12_scaled + bio10_15_scaled + 
               SnowMonths)^2 +
             #  SNDPPT # Not included, as the other two dictate the third
             
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = abundance,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

plot(a1)
simulationOutput_a1 <- simulateResiduals(fittedModel = a1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput_a1,quantreg = TRUE)
simulationOutput_a1b <- simulateResiduals(fittedModel = a1, refit = TRUE)
testOverdispersion(simulationOutput_a1b, alternative = "overdispersion", plot = TRUE)
# Not overdispersed
testZeroInflation(simulationOutput_a1, plot = TRUE, alternative = "more")
## But zeroinflated


abundance_model <- modelSimplification(model = a1, data = abundance, alpha = 0.05, optimizer = "bobyqa", Iters = 2e5)
save(abundance_model, file = file.path(models, "abundancemodel_full.rds"))



simulationOutput_a2 <- simulateResiduals(fittedModel = abundance_model, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput_a2,quantreg = TRUE)
simulationOutput_a2b <- simulateResiduals(fittedModel = abundance_model, refit = TRUE)
testOverdispersion(simulationOutput_a2, alternative = "overdispersion", plot = TRUE)
# Overdispersed
testZeroInflation(simulationOutput_a2, plot = TRUE, alternative = "more")
## But zeroinflated

