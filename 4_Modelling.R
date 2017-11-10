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

#sites <- read.csv(file.path(data_in, loadin))
sites <- read.csv("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\3_Data\\Sites_2017-11-09.csv")
rm(loadin)

#################################################
# 4. Set reference levels
#################################################

sites <- SiteLevels(sites) ## relevels all land use/habitat variables


#################################################
# 5. Investigating Climate Variables
#################################################

## RandomForest classification to determine variable importance
spR <- sites[!(is.na(sites$SpeciesRichness)),]
spR1 <- randomForest(SpeciesRichness ~ bio10_1 + bio10_5 + bio10_6 + bio10_12 + bio10_13 + bio10_14, data = spR, importance=TRUE)
# % Var explained: 55.35
plot(spR1)
varImpPlot(spR1, type = 2)
importance(spR1)

### Most important are bio10_14 nd then bio10_13

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

richness <- richness[complete.cases(richness$scalePH),]


table(richness$ESA)
richness_notinclude <- c("Needleleaf deciduous forest", "Tree open", "Herbaceous with spare tree/shrub ", 
                         "Sparse vegetation", "Cropland/Other vegetation mosaic", 
                         "Bare area (consolidated", "Paddy field", "Wetland/Herbaceous")

richness <- droplevels(richness[!(richness$ESA %in% richness_notinclude),])
richness$scalePH <- as.vector(scale(richness$ph_new))
richness$bio10_5_scaled <- scale(richness$bio10_5)
richness$bio10_13_scaled <- scale(richness$bio10_13)
richness$bio10_14_scaled <- scale(richness$bio10_14)

r1 <- glmer(SpeciesRichness ~  scalePH  + ESA + 
             scalePH:ESA + 
             (bio10_14_scaled + bio10_13_scaled + bio10_5_scaled)^2 + 
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = richness, family = poisson,
            control = glmerControl(optCtrl = list(maxfun = 100000), optimizer ="Nelder_Mead"))


summary(r1)

r2a <- update(r1, .~. -bio10_13_scaled:bio10_5_scaled)
anova(r1, r2a) # Not significnat
r2b <- update(r1, .~. -bio10_14_scaled:bio10_5_scaled)
anova(r1, r2b) # Not significant
r2c <- update(r1, .~. -bio10_14_scaled:bio10_13_scaled)
anova(r1, r2c) #Not significant

r2d <- update(r1, .~. -scalePH:ESA)
anova(r1, r2d) # Not significant


r3a <- update(r2d, .~. -bio10_13_scaled:bio10_5_scaled)
anova(r2d, r3a) # Not significnat
r3b <- update(r2d, .~. -bio10_14_scaled:bio10_5_scaled)
anova(r2d, r3b) # Not significant
r3c <- update(r2d, .~. -bio10_14_scaled:bio10_13_scaled)
anova(r2d, r3c) #Not significant

r4a <- update(r3a, .~. -bio10_14_scaled:bio10_5_scaled)
anova(r3a, r4a)#Not significant
r4b <- update(r3a, .~. -bio10_14_scaled:bio10_13_scaled)
anova(r3a, r4b)#Not significant


r5a <- update(r4a, .~. -bio10_14_scaled:bio10_13_scaled)
anova(r4a, r5a)#Not significant

r6a <- update(r5a, .~. -bio10_5_scaled)
anova(r5a, r6a) #Not significant
r6b <- update(r5a, .~. -bio10_13_scaled)
anova(r5a, r6b)  #Not significant
r6c <- update(r5a, .~. -bio10_14_scaled)
anova(r5a, r6c) #Not significant
r6d <- update(r5a, .~. -ESA)
anova(r5a, r6d) #Not significant
r6e <- update(r5a, .~. -scalePH)
anova(r5a, r6e) #Not significant


r7a <- update(r6c, .~. -bio10_5_scaled)
anova(r6c, r7a) #Not significant
r7b <- update(r6c, .~. -bio10_13_scaled)
anova(r6c, r7b) #Not significant
r7c <- update(r6c, .~. -ESA)
anova(r6c, r7c) #Not significant
r7d <- update(r6c, .~. -scalePH)
anova(r6c, r7d) #Not significant


r8a <- update(r7b, .~. -bio10_5_scaled)
anova(r7b, r8a) #Not significant
r8b <- update(r7b, .~. -ESA)
anova(r7b, r8b) #Not significant
r8c <- update(r7b, .~. -scalePH)
anova(r7b, r8b) #Not significant

r9a <- update(r8a, .~. -ESA)
anova(r8a, r9a) #Not significant
r9b <- update(r8a, .~. -scalePH)
anova(r8a, r9b) #Not significant

r10 <- update(r9a, .~. -scalePH)
anova(r9a, r10) ## Not significant

anova(r1, r10) ## Not sigificnat
## well that's a bit shit

#plot(sp_habitat)
#save(sp_habitat, file = file.path(models, "sp_habitat.rds"))

#################################################
# 5. Biomass
#################################################
biomass <- sites[complete.cases(sites$logBiomass),]
biomass <- droplevels(biomass[biomass$ESA != "Unknown",])
biomass$scalePH <- as.vector(scale(biomass$ph_new))


table(biomass$ESA)
biomass_notinclude <- c("Tree open", "Sparse vegetation", "Cropland/Other vegetation mosaic", "Urban", "Paddy field", "Herbaceous with spare tree/shrub")

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
abundance$scalePH <- scale(abundance$ph_new)

abundance_notinclude <- c("Needleleaf deciduous forest", "Sparse vegetation",
                          "Cropland/Other vegetation mosaic", "Bare area (consolidated", "Paddy field",
                          "Wetland/Herbaceous", "Herbaceous with spare tree/shrub")

abundance <- droplevels(abundance[!(abundance$ESA %in% abundance_notinclude),])
tapply(abundance$scalePH, abundance$ESA, summary)

abundance$bio10_5_scaled <- scale(abundance$bio10_5)
abundance$bio10_13_scaled <- scale(abundance$bio10_13)
abundance$bio10_14_scaled <- scale(abundance$bio10_14)

a1 <- lmer(logAbundance ~  scalePH  + ESA + 
             scalePH:ESA + 
             (bio10_14_scaled * bio10_13_scaled * bio10_5_scaled) + 
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = abundance)

plot(a1)


a2a <- update(a1, .~. -bio10_14_scaled:bio10_13_scaled:bio10_5_scaled)
anova(a1, a2a) # Not significant

a3a <- update(a2a, .~. -bio10_13_scaled:bio10_5_scaled)
anova(a2a, a3a) # Not significant
a3b <- update(a2a, .~. -bio10_14_scaled:bio10_5_scaled )
anova(a2a, a3b) # Not significant
a3c <- update(a2a, .~. -bio10_14_scaled:bio10_13_scaled)
anova(a2a, a3c) # Not significant

a4a <- update(a3c, .~. -bio10_13_scaled:bio10_5_scaled)
anova(a3c, a4a)# Not significant
a4b <- update(a3c, .~. -bio10_14_scaled:bio10_5_scaled)
anova(a3c, a4b)# Not significant

a5a <- update(a4b, .~. -bio10_13_scaled:bio10_5_scaled)
anova(a4b, a5a) # Significant

a6a <- update(a4b, .~. -scalePH:ESA)
anova(a4b, a6a) ## Significant

a7a <- update(a4b, .~. -bio10_14_scaled)
anova(a4b, a7a) ## Significant

summary(a4b)
####
## a1
####
save(a4b, file = file.path(models, "abundance_ESApH+CHELSA.rds.rds"))

