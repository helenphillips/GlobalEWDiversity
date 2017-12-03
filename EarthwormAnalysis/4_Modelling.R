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

## RandomForest classification to determine variable importance
## Species Richness
spR <- sites[!(is.na(sites$SpeciesRichness)),]
spR <- spR[!(is.na(spR$bio10_19)),]
spR1 <- randomForest(SpeciesRichness ~ bio10_1 + bio10_2 +bio10_3 +bio10_4 +bio10_5 
                     + bio10_6 +bio10_7 +bio10_8 +bio10_9 +bio10_10 +bio10_11
                     + bio10_12 + bio10_13 + bio10_14 + bio10_15 +bio10_16 +bio10_17 
                     +bio10_18 +bio10_19 , data = spR, importance=TRUE)
# % Var explained: 54.16
plot(spR1)
varImpPlot(spR1, type = 2)
importance(spR1)


## RandomForest classification to determine variable importance
## Abundance
Ab <- sites[!(is.na(sites$Site_Abundance)),]
Ab <- Ab[!(is.na(Ab$bio10_19)),]
Ab1 <- randomForest(Site_Abundance ~ bio10_1 + bio10_2 +bio10_3 +bio10_4 +bio10_5 
                     + bio10_6 +bio10_7 +bio10_8 +bio10_9 +bio10_10 +bio10_11
                     + bio10_12 + bio10_13 + bio10_14 + bio10_15 +bio10_16 +bio10_17 
                     +bio10_18 +bio10_19 , data = Ab, importance=TRUE)
# % Var explained: 36.35
plot(Ab1)
varImpPlot(Ab1, type = 2)
importance(Ab1)


## RandomForest classification to determine variable importance
## Biomass
Bio <- sites[!(is.na(sites$Site_WetBiomass)),]
Bio <- Bio[!(is.na(Bio$bio10_19)),]
Bio1 <- randomForest(Site_WetBiomass ~ bio10_1 + bio10_2 +bio10_3 +bio10_4 +bio10_5 
                    + bio10_6 +bio10_7 +bio10_8 +bio10_9 +bio10_10 +bio10_11
                    + bio10_12 + bio10_13 + bio10_14 + bio10_15 +bio10_16 +bio10_17 
                    +bio10_18 +bio10_19 , data = Bio, importance=TRUE)
# % Var explained: -2.43
plot(Bio1)
varImpPlot(Bio1, type = 2)
importance(Bio1)



## Checking for collinearity
plot(sites$bio10_14 ~ sites$bio10_13)
r1 <- lm(sites$bio10_14 ~ sites$bio10_13)
summary(r1) # low r1, with no random effect structure

##https://stats.stackexchange.com/questions/82984/how-to-test-and-avoid-multicollinearity-in-mixed-linear-model/142032
## That page said a cut of 4. We are at 0.3
x <- data.frame(sites$bio10_14,sites$bio10_13)
correl_dummy_df <- round(cor(x, use = "pair"), 2) 

## VIFs
source("MEE3_1_sm_Appendix_S1/HighstatLib.R")
z <- data.frame(sites$bio10_14,sites$bio10_13,sites$bio10_5)
corvif(z) ## All three seem fine

#################################################
# 4. Species Richness
#################################################

richness <- sites[complete.cases(sites$SpeciesRichness),]
richness <- droplevels(richness[richness$ESA != "Unknown",])
richness <- droplevels(richness[-which(richness$SpeciesRichness != round(richness$SpeciesRichness)),])

# richness <- richness[complete.cases(richness$scalePH),]

richness <- droplevels(richness[!(is.na(richness$PHIHOX)),])
richness <- droplevels(richness[!(is.na(richness$bio10_2)),]) ## 1781


table(richness$ESA)
richness_notinclude <- c("Needleleaf deciduous forest", 
                         "Sparse vegetation",  
                         "Bare area (consolidated", "Paddy field")

richness <- droplevels(richness[!(richness$ESA %in% richness_notinclude),]) ##  1764
summary(richness$PHIHOX)
richness$scalePH <- as.vector(scale(richness$PHIHOX))
richness$scaleCLYPPT <- scale(richness$CLYPPT)
richness$scaleSLTPPT <- scale(richness$SLTPPT)
richness$scaleCECSOL <- scale(richness$CECSOL)
richness$scaleORCDRC <- scale(richness$ORCDRC)


# Three most important climate variables
richness$bio10_2_scaled <- scale(richness$bio10_2)
richness$bio10_15_scaled <- scale(richness$bio10_15)
richness$bio10_5_scaled <- scale(richness$bio10_5)

r1 <- glmer(SpeciesRichness ~  ESA + (scalePH  + 
             scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_2_scaled + bio10_15_scaled + bio10_5_scaled)^2 + 
             #  SNDPPT # Not included, as the other two dictate the third
              
              # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
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


#################################################
# 5. Biomass
#################################################
biomass <- sites[complete.cases(sites$logBiomass),]
biomass <- droplevels(biomass[biomass$ESA != "Unknown",])

biomass <- droplevels(biomass[!(is.na(biomass$PHIHOX)),])
biomass <- droplevels(biomass[!(is.na(biomass$bio10_2)),]) ## 1011


table(biomass$ESA)
biomass_notinclude <- c("Tree open", "Sparse vegetation", 
                        "Urban", "Paddy field", "Herbaceous with spare tree/shrub")

biomass <- droplevels(biomass[!(biomass$ESA %in% biomass_notinclude),]) ##  969
summary(biomass$PHIHOX)
biomass$scalePH <- as.vector(scale(biomass$PHIHOX))
biomass$scaleCLYPPT <- scale(biomass$CLYPPT)
biomass$scaleSLTPPT <- scale(biomass$SLTPPT)
biomass$scaleCECSOL <- scale(biomass$CECSOL)
biomass$scaleORCDRC <- scale(biomass$ORCDRC)



biomass <- droplevels(biomass[!(biomass$ESA %in% biomass_notinclude),])

biomass$bio10_5_scaled <- scale(biomass$bio10_5)
biomass$bio10_13_scaled <- scale(biomass$bio10_13)
biomass$bio10_14_scaled <- scale(biomass$bio10_14)



b1 <- lmer(logBiomass ~  scalePH  + ESA + 
              scalePH:ESA + 
             (bio10_14_scaled * bio10_13_scaled * bio10_5_scaled) + 
              # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
            #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
            # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
            (1|file/Study_Name), data = biomass)

b2a <- update(b1, .~. -bio10_14_scaled:bio10_13_scaled:bio10_5_scaled)
anova(b1, b2a) ## Not significant

b3a <- update(b2a, .~. -bio10_13_scaled:bio10_5_scaled)
anova(b2a, b3a) ## Not Significant

b3b <- update(b2a, .~. -bio10_14_scaled:bio10_5_scaled)
anova(b2a, b3b) ## Not Significant

b3c <- update(b2a, .~. -bio10_14_scaled:bio10_13_scaled)
anova(b2a, b3c) ## Significant


b4a <- update(b3b, .~. -bio10_13_scaled:bio10_5_scaled)
anova(b3b, b4a) ## Not Significant

b4b <- update(b3b, .~. -bio10_14_scaled:bio10_13_scaled)
anova(b3b, b4b) ## Significant


b5a <- update(b4a, .~. -bio10_14_scaled:bio10_13_scaled)
anova(b4a, b5a) ##  Significant

summary(b4a)


b6a <- update(b4a, .~. -scalePH:ESA)
anova(b5a, b6a) ## Significant


b7a <- update(b4a, .~. - bio10_5_scaled)
anova(b4a, b7a) ## Significant

summary(b4a)

######
## b2a
######
save(b4a, file = file.path(models, "biomass_ESApH+CHELSA.rds"))


#################################################
# 6. Abundance
#################################################
hist(sites$Site_Abundance)
hist(sites$logAbundance)
abundance <- sites[complete.cases(sites$Site_Abundance),]
abundance <- droplevels(abundance[abundance$ESA != "Unknown",])

abundance <- droplevels(abundance[!(is.na(abundance$PHIHOX)),])
abundance <- droplevels(abundance[!(is.na(abundance$bio10_2)),]) ## 2091


table(abundance$ESA)
abundance_notinclude <- c("Needleleaf deciduous forest", "Sparse vegetation",
                         "Bare area (consolidated", "Paddy field")

abundance <- droplevels(abundance[!(abundance$ESA %in% abundance_notinclude),])
tapply(abundance$scalePH, abundance$ESA, summary)


abundance$scalePH <- as.vector(scale(abundance$PHIHOX))
abundance$scaleCLYPPT <- scale(abundance$CLYPPT)
abundance$scaleSLTPPT <- scale(abundance$SLTPPT)
abundance$scaleCECSOL <- scale(abundance$CECSOL)
abundance$scaleORCDRC <- scale(abundance$ORCDRC)

abundance$bio10_2_scaled <- scale(abundance$bio10_2)
abundance$bio10_10_scaled <- scale(abundance$bio10_10)
abundance$bio10_18_scaled <- scale(abundance$bio10_18)

a1 <- lmer(logAbundance ~  ESA + (scalePH  + 
                                    scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_2_scaled + bio10_10_scaled + bio10_18_scaled)^2 + 
             #  SNDPPT # Not included, as the other two dictate the third
             
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = abundance, family = "gaussian"
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

plot(a1)
abundance_model <- modelSimplification(model = a1, data = abundance, alpha = 0.05, optimizer = "bobyqa", Iters = 2e5)
save(abundance_model, file = file.path(models, "abundancemodel_full.rds"))
