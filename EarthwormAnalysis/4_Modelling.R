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
library(randomForest)
source("Functions/FormatData.R")
source("Functions/lme4_ModellingFunctions.R")
source("Functions/ModelSimplification.R")

#################################################
# 2. Loading in variables
#################################################

data_in <-"3_Data"

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
x <- data.frame(sites$bio10_1,sites$bio10_4,sites$bio10_7,sites$bio10_12,sites$bio10_15)
correl_dummy_df <- round(cor(x, use = "pair"), 2) 

## VIFs
source("MEE3_1_sm_Appendix_S1/HighstatLib.R")
corvif(x) ## There might be an issue here
y <- data.frame(sites$bio10_1,sites$bio10_4,sites$bio10_12,sites$bio10_15)
corvif(y)

#################################################
# 4. Species Richness
#################################################

richness <- sites[complete.cases(sites$SpeciesRichness),] #3119
richness <- droplevels(richness[richness$ESA != "Unknown",])
richness <- droplevels(richness[-which(richness$SpeciesRichness != round(richness$SpeciesRichness)),])

# richness <- richness[complete.cases(richness$scalePH),]

richness <- droplevels(richness[!(is.na(richness$PHIHOX)),])
richness <- droplevels(richness[!(is.na(richness$bio10_2)),]) ## 2645


table(richness$ESA)
richness_notinclude <- c("Needleleaf deciduous forest", 
                         "Sparse vegetation",  "Cropland/Other vegetation mosaic",
                         "Bare area (consolidated", "Paddy field", "Wetland/Herbaceous")

richness <- droplevels(richness[!(richness$ESA %in% richness_notinclude),]) ##  2623
summary(richness$PHIHOX)
richness$scalePH <- as.vector(scale(richness$PHIHOX))
richness$scaleCLYPPT <- scale(richness$CLYPPT)
richness$scaleSLTPPT <- scale(richness$SLTPPT)
richness$scaleCECSOL <- scale(richness$CECSOL)
richness$scaleORCDRC <- scale(richness$ORCDRC)


# Three most important climate variables
richness$bio10_1_scaled <- scale(richness$bio10_1)
richness$bio10_4_scaled <- scale(richness$bio10_4)
# richness$bio10_7_scaled <- scale(richness$bio10_7)
richness$bio10_12_scaled <- scale(richness$bio10_12)
richness$bio10_15_scaled <- scale(richness$bio10_15)
## Save the data
write.csv(richness, file = file.path(data_out, paste("sitesRichness_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

r1 <- glmer(SpeciesRichness ~  ESA + (scalePH  + 
             scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_1_scaled + bio10_4_scaled + 
                bio10_12_scaled + bio10_15_scaled)^2 + 
             #  SNDPPT # Not included, as the other two dictate the third
              
             #  (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = richness, family = poisson,
            control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))


summary(r1)



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

#################################################
# 5. Biomass
#################################################
biomass <- sites[complete.cases(sites$logBiomass),]
biomass <- droplevels(biomass[biomass$ESA != "Unknown",])

biomass <- droplevels(biomass[!(is.na(biomass$PHIHOX)),])
biomass <- droplevels(biomass[!(is.na(biomass$bio10_2)),]) ## 1605


table(biomass$ESA)
biomass_notinclude <- c("Tree open", "Sparse vegetation", "Cropland/Other vegetation mosaic",
                        "Urban", "Paddy field")

biomass <- droplevels(biomass[!(biomass$ESA %in% biomass_notinclude),]) ##  1585
summary(biomass$PHIHOX)
biomass$scalePH <- as.vector(scale(biomass$PHIHOX))
biomass$scaleCLYPPT <- scale(biomass$CLYPPT)
biomass$scaleSLTPPT <- scale(biomass$SLTPPT)
biomass$scaleCECSOL <- scale(biomass$CECSOL)
biomass$scaleORCDRC <- scale(biomass$ORCDRC)

biomass$bio10_1_scaled <- scale(biomass$bio10_1)
biomass$bio10_4_scaled <- scale(biomass$bio10_4)
biomass$bio10_7_scaled <- scale(biomass$bio10_7)
biomass$bio10_12_scaled <- scale(biomass$bio10_12)
biomass$bio10_15_scaled <- scale(biomass$bio10_15)

## Save the data
write.csv(biomass, file = file.path(data_out, paste("sitesBiomass_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)



b1 <- lmer(logBiomass ~  ESA + (scalePH  + 
            scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_1_scaled + bio10_4_scaled + 
                bio10_12_scaled + bio10_15_scaled)^2 +
             #  SNDPPT # Not included, as the other two dictate the third
             
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = biomass,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))
plot(b1)

biomass_model <- modelSimplification(model = b1, data = biomass, alpha = 0.05, optimizer = "bobyqa", Iters = 2e5)
save(biomass_model, file = file.path(models, "biomassmodel_full.rds"))
# load(file.path(models, "richnessmodel_full.rds"))


#################################################
# 6. Abundance
#################################################
hist(sites$Site_Abundance)
hist(sites$logAbundance)
abundance <- sites[complete.cases(sites$Site_Abundance),]
abundance <- droplevels(abundance[abundance$ESA != "Unknown",]) #3455

abundance <- droplevels(abundance[!(is.na(abundance$PHIHOX)),])
abundance <- droplevels(abundance[!(is.na(abundance$bio10_2)),]) ##  3329


table(abundance$ESA)
abundance_notinclude <- c("Needleleaf deciduous forest", "Sparse vegetation", "Cropland/Other vegetation mosaic", 
                         "Bare area (consolidated", "Bare area (unconsolidated",  "Paddy field", "Wetland/Herbaceous")

abundance <- droplevels(abundance[!(abundance$ESA %in% abundance_notinclude),]) #  3284
tapply(abundance$scalePH, abundance$ESA, summary)


abundance$scalePH <- as.vector(scale(abundance$PHIHOX))
abundance$scaleCLYPPT <- scale(abundance$CLYPPT)
abundance$scaleSLTPPT <- scale(abundance$SLTPPT)
abundance$scaleCECSOL <- scale(abundance$CECSOL)
abundance$scaleORCDRC <- scale(abundance$ORCDRC)

abundance$bio10_1_scaled <- scale(abundance$bio10_1)
abundance$bio10_4_scaled <- scale(abundance$bio10_4)
abundance$bio10_7_scaled <- scale(abundance$bio10_7)
abundance$bio10_12_scaled <- scale(abundance$bio10_12)
abundance$bio10_15_scaled <- scale(abundance$bio10_15)


## Save the data
write.csv(abundance, file = file.path(data_out, paste("sitesAbundance_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)


a1 <- lmer(logAbundance ~  ESA + (scalePH  + 
                                    scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_1_scaled + bio10_4_scaled + 
                bio10_12_scaled + bio10_15_scaled)^2 +
             #  SNDPPT # Not included, as the other two dictate the third
             
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = abundance,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

plot(a1)
abundance_model <- modelSimplification(model = a1, data = abundance, alpha = 0.05, optimizer = "bobyqa", Iters = 2e5)
save(abundance_model, file = file.path(models, "abundancemodel_full.rds"))
