########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}


if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\Users\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
  
}



#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
library(lme4)
library(glmmTMB)
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

data_in <-"16_Data"

files <- list.files(file.path(data_in))
files <- files[grep("Sites", files)]
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

if(!dir.exists("18_Data")){
  dir.create("18_Data")
}

data_out <- "18_Data"

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


## This data has already been subset to the data that was used
## in the original manuscript


#################################################
# 4. Set reference levels
#################################################

sites <- SiteLevels(sites) ## relevels all land use/habitat variables


#################################################
# 5. Reorder data frame, column with FG
#################################################
## Remove functional group diversity measures to use as idvars
idvars <- names(sites)[-(114:ncol(sites))]

m_sites <- melt(sites, id.vars = idvars)

#################################################
# 6. Split into three diversity measures
#################################################
biomass <- m_sites[grep("biomass", m_sites$variable),]
biomass <- droplevels(biomass[which(!(is.na(biomass$value))),])
nobiomass <- aggregate(biomass$value ~ biomass$Study_Name, FUN = var)
nobiomass <- nobiomass[which(nobiomass[,2] == 0 | is.na(nobiomass[,2])),1]
biomass <- droplevels(biomass[!(biomass$Study_Name %in% nobiomass),]) # 8817 # 8035


biomass <- droplevels(biomass[biomass$ESA != "Unknown",]) # 7927 # 7960

#  biomass <- biomass[complete.cases(biomass$value),]

biomass <- droplevels(biomass[!(is.na(biomass$bio10_15)),]) ## 3365
biomass <- droplevels(biomass[!(is.na(biomass$OCFinal)),]) ## 3364 # 7960
biomass <- droplevels(biomass[!(is.na(biomass$phFinal)),]) ## 3364
biomass <- droplevels(biomass[!(is.na(biomass$SnowMonths_cat)),]) ##  3361
biomass <- droplevels(biomass[!(is.na(biomass$Aridity)),]) ##  3357 # 7955



table(biomass$ESA, biomass$variable)

biomass_notinclude <- c("Cropland/Other vegetation mosaic")
biomass <- droplevels(biomass[!(biomass$ESA %in% biomass_notinclude),]) ##   7772

fg_biomass <- biomass

# Scale variables
biomass <- scaleVariables(biomass)

## Save the data
write.csv(biomass, file = file.path(data_out, paste("sitesFGBiomass_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

## Check the data

randomIDs <- sample(biomass$ID, size = 10)
testB <- droplevels(biomass[biomass$ID %in% randomIDs,])
# Manual check of original data

# findVariables(biomass, VIFThreshold = 3)
ind <- df_variables(biomass)
dat <- biomass[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# Remove
# Bio 7
# Bio 1
# Aridity
# Petyr


############################### Abundance
abundance <- m_sites[grep("abundance", m_sites$variable),]
abundance <- droplevels(abundance[which(!(is.na(abundance$value))),])
noabundance <- aggregate(abundance$value ~ abundance$Study_Name, FUN = var)
noabundance <- noabundance[which(noabundance[,2] == 0 | is.na(noabundance[,2])),1]
abundance <- droplevels(abundance[!(abundance$Study_Name %in% noabundance),]) # 25384 #  25995

abundance <- droplevels(abundance[abundance$ESA != "Unknown",]) #6759 # 25065


# abundance <- droplevels(abundance[!(is.na(abundance$PHIHOX)),])
abundance <- droplevels(abundance[!(is.na(abundance$bio10_15)),]) ##   
abundance <- droplevels(abundance[!(is.na(abundance$OCFinal)),]) ##  
abundance <- droplevels(abundance[!(is.na(abundance$phFinal)),]) ##  6731
abundance <- droplevels(abundance[!(is.na(abundance$SnowMonths_cat)),]) ##  6657
abundance <- droplevels(abundance[!(is.na(abundance$Aridity)),]) ##  6576 # 24590


table(abundance$ESA, abundance$variable)
abundance_notinclude <- c("Needleleaf deciduous forest", "Tree open", "Sparse vegetation", "Cropland/Other vegetation mosaic",
                         "Wetland/Herbaceous", "Urban",
                          "Water bodies")

abundance <- droplevels(abundance[!(abundance$ESA %in% abundance_notinclude),]) #  25107 # 24490

fg_abundance <- abundance

# Scale variables
abundance <- scaleVariables(abundance)



write.csv(abundance, file = file.path(data_out, paste("sitesFGAbundance_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

# findVariables(biomass, VIFThreshold = 3)
ind <- df_variables(abundance)
dat <- abundance[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# Remove bio 4
# Remove 1
# Remove 12
# Remove pet_sd
## All ok



#################################################
# 7. Biomass Modelling
#################################################


# biomass$value[which(biomass$value < 0)] <- 0
# For now

biomass$logValue <- log(biomass$value + 1)
hist(biomass$logValue)

b1 <- glmmTMB(logBiomass ~  (ESA * variable) + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                (bio10_7_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePET + SnowMonths_cat)^2 + 
                scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
                scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                (1|file/Study_Name), data = biomass, 
              ziformula = ~ESA + variable + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL +
                bio10_7_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePET + SnowMonths_cat,
              control = glmmTMBControl(optCtrl = list(iter.max = 2e5,eval.max=2e5))) #, optimizer ="bobyqa"))

simulationOutput_bm <- simulateResiduals(fittedModel = b1, n = 250)
plot(simulationOutput_bm,quantreg = TRUE)
# so ignore the res versus pred plot, but look at others
testDispersion(simulationOutput_bm, alternative = "greater", plot = TRUE)
# A bit overdispersed
testZeroInflation(simulationOutput_bm, plot = TRUE, alternative = "greater")
## But zeroinflated


# b1 <- lmer(logValue ~  (ESA * ) + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
#              (bio10_4_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePETSD + SnowMonths_cat)^2 + 
#              scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
#              scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
#              ScalePETSD:bio10_12_scaled + ScalePETSD:bio10_15_scaled +
#              (1|file/Study_Name), data = biomass,
#            control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))
# 

biomass_model <- modelSimplificationAIC(model = b1, data = biomass, optimizer = "bobyqa", Iters = 2e5)
save(biomass_model, file = file.path(models, "biomassmodel_functionalgroups_revised.rds"))
# load(file.path(models, "biomassmodel_functionalgroups.rds"))

############### abundance
abundance$logValue <- log(abundance$value + 1)
hist(abundance$logValue)

a1 <- lmer(logValue ~  (ESA * variable) + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                ScalePET)^2 +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
             scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET + 
             scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity + 
             (1|file/Study_Name), data = abundance,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

abundance_model <- modelSimplificationAIC(model = a1, data = abundance, optimizer = "bobyqa", Iters = 2e5)
save(abundance_model, file = file.path(models, "abundancemodel_functionalgroups_revised.rds"))
# load(file.path(models, "abundancemodel_functionalgroups.rds"))


########################
# NOTHING PAST THIS POINT HAS BEEN UPDATED!!
########################


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



