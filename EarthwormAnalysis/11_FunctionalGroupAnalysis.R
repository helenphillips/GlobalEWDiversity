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
library(reshape)
source("Functions/FormatData.R")
source("Functions/CorvifVariablePicker.R")
source("Functions/lme4_ModellingFunctions.R")
source("Functions/ModelSimplification.R")
source("MEE3_1_sm_Appendix_S1/HighstatLib.R")

#################################################
# 2. Loading in variables
#################################################

data_in <-"10_Data"

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

if(!dir.exists("11_Data")){
  dir.create("11_Data")
}

data_out <- "11_Data"

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
# 5. Reorder data frame, column with FG
#################################################

idvars <- names(sites)[-(111:ncol(sites))]

m_sites <- melt(sites, id.vars = idvars)

#################################################
# 6. Split into three diversity measures
#################################################
biomass <- m_sites[grep("biomass", m_sites$variable),]
biomass <- droplevels(biomass[which(!(is.na(biomass$value))),])
nobiomass <- aggregate(biomass$value ~ biomass$Study_Name, FUN = var)
nobiomass <- nobiomass[which(nobiomass[,2] == 0 | is.na(nobiomass[,2])),1]
biomass <- droplevels(biomass[!(biomass$Study_Name %in% nobiomass),]) # 8817


biomass <- droplevels(biomass[biomass$ESA != "Unknown",]) # 7927

#  biomass <- biomass[complete.cases(biomass$value),]

biomass <- droplevels(biomass[!(is.na(biomass$bio10_15)),]) ## 3365
biomass <- droplevels(biomass[!(is.na(biomass$OCFinal)),]) ## 3364
biomass <- droplevels(biomass[!(is.na(biomass$phFinal)),]) ## 3364
biomass <- droplevels(biomass[!(is.na(biomass$SnowMonths_cat)),]) ##  3361
biomass <- droplevels(biomass[!(is.na(biomass$Aridity)),]) ##  3357



table(biomass$ESA, biomass$variable)

biomass_notinclude <- c("Cropland/Other vegetation mosaic")
biomass <- droplevels(biomass[!(biomass$ESA %in% biomass_notinclude),]) ##   7772

fg_biomass <- biomass

# Scale variables
biomass <- scaleVariables(biomass)

## Save the data
write.csv(biomass, file = file.path(data_out, paste("sitesFGBiomass_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

# findVariables(biomass, VIFThreshold = 3)
ind <- df_variables(biomass)
dat <- biomass[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# Remove bio 7
# REmove 1
# Remove petyr
# Remove Aridity
# Ok

############################### Abundance
abundance <- m_sites[grep("abundance", m_sites$variable),]
abundance <- droplevels(abundance[which(!(is.na(abundance$value))),])
noabundance <- aggregate(abundance$value ~ abundance$Study_Name, FUN = var)
noabundance <- noabundance[which(noabundance[,2] == 0 | is.na(noabundance[,2])),1]
abundance <- droplevels(abundance[!(abundance$Study_Name %in% noabundance),]) # 25384

abundance <- droplevels(abundance[abundance$ESA != "Unknown",]) #6759


# abundance <- droplevels(abundance[!(is.na(abundance$PHIHOX)),])
abundance <- droplevels(abundance[!(is.na(abundance$bio10_15)),]) ##   
abundance <- droplevels(abundance[!(is.na(abundance$OCFinal)),]) ##  
abundance <- droplevels(abundance[!(is.na(abundance$phFinal)),]) ##  6731
abundance <- droplevels(abundance[!(is.na(abundance$SnowMonths_cat)),]) ##  6657
abundance <- droplevels(abundance[!(is.na(abundance$Aridity)),]) ##  6576


table(abundance$ESA, abundance$variable)
abundance_notinclude <- c("Needleleaf deciduous forest", "Tree open", "Sparse vegetation", "Cropland/Other vegetation mosaic",
                         "Wetland/Herbaceous",
                          "Water bodies")

abundance <- droplevels(abundance[!(abundance$ESA %in% abundance_notinclude),]) #  25107

fg_abundance <- abundance

# Scale variables
abundance <- scaleVariables(abundance)



write.csv(abundance, file = file.path(data_out, paste("sitesFGAbundance_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

# findVariables(biomass, VIFThreshold = 3)
ind <- df_variables(abundance)
dat <- abundance[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# Remove bio 7
# Remove 1
# Remove 12
# Remove pet_sd
## All ok



#################################################
# 7. Biomass Modelling
#################################################


biomass$value[which(biomass$value < 0)] <- 0
# For now

biomass$logValue <- log(biomass$value + 1)
hist(biomass$logValue)

b1 <- lmer(logValue ~  (ESA * variable) + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
             (bio10_12_scaled  + bio10_15_scaled + SnowMonths_cat)^2 + 
             scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
             #  SNDPPT # Not included, as the other two dictate the third
             
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = biomass,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

biomass_model <- modelSimplificationAIC(model = b1, data = biomass, optimizer = "bobyqa", Iters = 2e5)
save(biomass_model, file = file.path(models, "biomassmodel_functionalgroups.rds"))
# load(file.path(models, "biomassmodel_functionalgroups.rds"))

############### abundance
abundance$logValue <- log(abundance$value + 1)
hist(abundance$logValue)


a1 <- lmer(logValue ~  (ESA * variable) + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
             (bio10_4_scaled  + bio10_15_scaled + SnowMonths_cat +scaleAridity + ScalePET)^2 + 
             scaleCLYPPT:bio10_4_scaled + scaleSLTPPT:bio10_4_scaled +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
             scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
             scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
             #  SNDPPT # Not included, as the other two dictate the third
             
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = abundance,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

abundance_model <- modelSimplificationAIC(model = a1, data = abundance, optimizer = "bobyqa", Iters = 2e5)
save(abundance_model, file = file.path(models, "abundancemodel_functionalgroups.rds"))
# load(file.path(models, "abundancemodel_functionalgroups.rds"))

############## Richness
richness <- m_sites[grep("richness", m_sites$variable),]
richness <- droplevels(richness[which(!(is.na(richness$value))),])
norichness <- aggregate(richness$value ~ richness$Study_Name, FUN = mean) # often no variation, so use mean
norichness <- norichness[which(norichness[,2] == 0 | is.na(norichness[,2])),1]
richness <- droplevels(richness[!(richness$Study_Name %in% norichness),]) # 9481


richness <- droplevels(richness[richness$ESA != "Unknown",]) # 8749

#  biomass <- biomass[complete.cases(biomass$value),]

richness <- droplevels(richness[!(is.na(richness$bio10_15)),]) ## 
richness <- droplevels(richness[!(is.na(richness$OCFinal)),]) ## 
richness <- droplevels(richness[!(is.na(richness$phFinal)),]) ## 
richness <- droplevels(richness[!(is.na(richness$SnowMonths_cat)),]) ##  
richness <- droplevels(richness[!(is.na(richness$Aridity)),]) ##  8505



table(richness$ESA, richness$variable)

richness_notinclude <- c("Needleleaf deciduous forest", "Tree open", "Sparse vegetation",
                        "Cropland/Other vegetation mosaic", "Bare area (consolidated", "Wetland/Herbaceous", "Water bodies")
richness <- droplevels(richness[!(richness$ESA %in% richness_notinclude),]) ##   8426

fg_richness <- richness

# Scale variables
richness <- scaleVariables(richness)


## Save the data
write.csv(richness, file = file.path(data_out, paste("sitesFGRichness_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)


# findVariables(biomass, VIFThreshold = 3)
ind <- df_variables(richness)
dat <- richness[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# Remove 7
# Remove 1
# Remove 12
# Remove petsd
# All ok

##### Modelling
## This is now done in a separate script on the cluster

# r1 <- glmer(value ~  (ESA * variable) + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
#              (bio10_4_scaled  + bio10_15_scaled + SnowMonths_cat +scaleAridity + ScalePET)^2 + 
#              scaleCLYPPT:bio10_4_scaled + scaleSLTPPT:bio10_4_scaled +
#              scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
#              scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
#              scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
#              (1|file/Study_Name), data = richness, family="poisson",
#            control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))
# 
# richness_model <- modelSimplificationAIC(model = r1, data = richness, optimizer = "bobyqa", Iters = 2e5)
# save(richness_model, file = file.path(models, "richnessmodel_functionalgroups.rds"))
# load(file.path(models, "abundancemodel_functionalgroups.rds"))


#################################################
# SPLIT BY FUNCTIONAL GROUP
#################################################

epi_biomass <- fg_biomass[grep("Epi", fg_biomass$variable),]
epi_biomass <- scaleVariables(epi_biomass)

# findVariables(biomass, VIFThreshold = 3)
ind <- df_variables(epi_biomass)
dat <- epi_biomass[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)
# Remove 7
# Remove 1
# Remove pet
# Remove aridity
# All ok


epi_biomass$value[which(epi_biomass$value < 0)] <- 0
# For now

epi_biomass$logValue <- log(epi_biomass$value + 1)
hist(epi_biomass$logValue)

epi_b1 <- lmer(logValue ~  ESA + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
             (bio10_4_scaled + bio10_12_scaled  + bio10_15_scaled + SnowMonths_cat + ScalePETSD)^2 + 
             scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
            scaleCLYPPT:ScalePETSD + scaleSLTPPT:ScalePETSD +
             #  SNDPPT # Not included, as the other two dictate the third
             
             # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
             
             # HabitatCover + 
             #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
             # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
             (1|file/Study_Name), data = epi_biomass,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

epi_biomass_model <- modelSimplificationAIC(model = epi_b1, data = epi_biomass, optimizer = "bobyqa", Iters = 2e5)
save(epi_biomass_model, file = file.path(models, "biomassmodel_epifunctionalgroups.rds"))


##############
endo_biomass <- fg_biomass[grep("Endo", fg_biomass$variable),]
endo_biomass <- scaleVariables(endo_biomass) # Function from FormatData.R
  
vars <- endo_biomass[,c(grep("scale", names(endo_biomass), ignore.case = TRUE))]
v <- findVariables(df = vars, VIFThreshold = 3)

# [1] "ScalePETSD"      "scalePH"         "scaleCLYPPT"     "scaleSLTPPT"    
# [5] "scaleCECSOL"     "scaleORCDRC"     "bio10_4_scaled"  "bio10_12_scaled"
# [9] "bio10_15_scaled"
endo_biomass$logValue <- log(endo_biomass$value + 1)

endo_b1 <- lmer(logValue ~  ESA + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                 (bio10_4_scaled + bio10_12_scaled  + bio10_15_scaled + SnowMonths_cat + ScalePETSD)^2 + 
                 scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
                 scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                 scaleCLYPPT:ScalePETSD + scaleSLTPPT:ScalePETSD +
                 #  SNDPPT # Not included, as the other two dictate the third
                 
                 # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
                 
                 # HabitatCover + 
                 #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
                 # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
                 (1|file/Study_Name), data = endo_biomass,
               control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

endo_biomass_model <- modelSimplificationAIC(model = endo_b1, data = endo_biomass, optimizer = "bobyqa", Iters = 2e5)
save(endo_biomass_model, file = file.path(models, "biomassmodel_endofunctionalgroups.rds"))


############
ane_biomass <- fg_biomass[grep("Ane", fg_biomass$variable),]
vars <- ane_biomass[,df_variables(ane_biomass)]

ane_biomass <- scaleVariables(ane_biomass)

v <- findVariables(df = vars, VIFThreshold = 3)

# [1] "ScalePETSD"      "scalePH"         "scaleCLYPPT"     "scaleSLTPPT"    
# [5] "scaleCECSOL"     "scaleORCDRC"     "bio10_4_scaled"  "bio10_12_scaled"
# [9] "bio10_15_scaled"
ane_biomass$logValue <- log(ane_biomass$value + 1)

ane_b1 <- lmer(logValue ~  ESA + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                  (bio10_4_scaled + bio10_12_scaled  + bio10_15_scaled + SnowMonths_cat + ScalePETSD)^2 + 
                  scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
                  scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                  scaleCLYPPT:ScalePETSD + scaleSLTPPT:ScalePETSD +
                  #  SNDPPT # Not included, as the other two dictate the third
                  
                  # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
                  
                  # HabitatCover + 
                  #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
                  # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
                  (1|file/Study_Name), data = ane_biomass,
                control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

ane_biomass_model <- modelSimplificationAIC(model = ane_b1, data = ane_biomass, optimizer = "bobyqa", Iters = 2e5)
save(ane_biomass_model, file = file.path(models, "biomassmodel_anefunctionalgroups.rds"))

########################################################3
## ABUNDANCE
#######################################################
#epi 
epi_abundance <- fg_abundance[grep("Epi", fg_abundance$variable),]
vars <- epi_abundance[,df_variables(epi_abundance)]

epi_abundance <- scaleVariables(epi_abundance)

v <- findVariables(df = vars, VIFThreshold = 3)


epi_abundance$logValue <- log(epi_abundance$value + 1)

epi_a1 <- lmer(logValue ~  ESA + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                  (bio10_4_scaled + bio10_15_scaled + SnowMonths_cat + ScalePET + scaleAridity)^2 + 
                  scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                  scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
                 scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
                  #  SNDPPT # Not included, as the other two dictate the third
                  
                  # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
                  
                  # HabitatCover + 
                  #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
                  # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
                  (1|file/Study_Name), data = epi_abundance,
                control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

epi_abundance_model <- modelSimplificationAIC(model = epi_a1, data = epi_abundance, optimizer = "bobyqa", Iters = 2e5)
save(epi_abundance_model, file = file.path(models, "abundancemodel_epifunctionalgroups.rds"))

## Endo
endo_abundance <- fg_abundance[grep("Endo", fg_abundance$variable),]
vars <- endo_abundance[,df_variables(endo_abundance)]

endo_abundance <- scaleVariables(endo_abundance)

v <- findVariables(df = vars, VIFThreshold = 3)


endo_abundance$logValue <- log(endo_abundance$value + 1)

endo_a1 <- lmer(logValue ~  ESA + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                  (bio10_4_scaled + bio10_15_scaled + SnowMonths_cat + ScalePET + scaleAridity)^2 + 
                  scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
                  scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                  scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
                  #  SNDPPT # Not included, as the other two dictate the third
                  
                  # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
                  
                  # HabitatCover + 
                  #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
                  # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
                  (1|file/Study_Name), data = endo_abundance,
                control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

endo_abundance_model <- modelSimplificationAIC(model = endo_a1, data = endo_abundance, optimizer = "bobyqa", Iters = 2e5)
save(endo_abundance_model, file = file.path(models, "abundancemodel_endofunctionalgroups.rds"))

## ane
ane_abundance <- fg_abundance[grep("Ane", fg_abundance$variable),]
vars <- ane_abundance[,df_variables(ane_abundance)]

ane_abundance <- scaleVariables(ane_abundance)

v <- findVariables(df = vars, VIFThreshold = 3)


ane_abundance$logValue <- log(ane_abundance$value + 1)

ane_a1 <- lmer(logValue ~  ESA + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                  (bio10_4_scaled + bio10_15_scaled + SnowMonths_cat + ScalePET + scaleAridity)^2 + 
                  scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
                  scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                  scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
                  #  SNDPPT # Not included, as the other two dictate the third
                  
                  # (Latitude__decimal_degrees * Longitude__Decimal_Degrees) +
                  
                  # HabitatCover + 
                  #   Soil_Organic_Matter__percent + # Organic_Carbon__percent +
                  # ph_new:HabitatCover + Organic_Carbon__percent:HabitatCover +
                  (1|file/Study_Name), data = ane_abundance,
                control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

ane_abundance_model <- modelSimplificationAIC(model = ane_a1, data = ane_abundance, optimizer = "bobyqa", Iters = 2e5)
save(ane_abundance_model, file = file.path(models, "abundancemodel_anefunctionalgroups.rds"))


########################################################3
## RICHNESS
#######################################################

# This is all done in a separate scrips on the cluster