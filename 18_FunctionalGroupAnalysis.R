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

m_sites <- melt(sites, id.vars = idvars) # 153856
## lots of rows, but its the three metrics, for five diff functional groups, but functional richness

#################################################
# 6. Split into three diversity measures
#################################################
biomass <- m_sites[grep("biomass", m_sites$variable),]
biomass <- droplevels(biomass[which(!(is.na(biomass$value))),]) # 28395
# 28395 / 5 # 5679



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
testB <-testB[,c("newID", "variable", "value")]


# findVariables(biomass, VIFThreshold = 3)
ind <- df_variables(biomass)
dat <- biomass[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# Remove
# Bio 4
# Bio 1
# Aridity
# Petyr


############################### Abundance
abundance <- m_sites[grep("abundance", m_sites$variable),] # 48080
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

# Remove petyr
# remove aridity
## All ok



#################################################
# 7. Biomass Modelling
#################################################

biomass$logValue <- log(biomass$value + 1)
hist(biomass$logValue)

b1 <- glmmTMB(logValue ~  (ESA * variable) + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                (bio10_7_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePETSD + SnowMonths_cat)^2 + 
                scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
                scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                (1|file/Study_Name), data = biomass, 
              ziformula = ~ESA + variable + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL +
                bio10_7_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePETSD + SnowMonths_cat,
              control = glmmTMBControl(optCtrl = list(iter.max = 2e5,eval.max=2e5))) #, optimizer ="bobyqa"))
beep()
simulationOutput_bm <- simulateResiduals(fittedModel = b1, n = 250)
plot(simulationOutput_bm,quantreg = TRUE)
# so ignore the res versus pred plot, but look at others
testDispersion(simulationOutput_bm, alternative = "greater", plot = TRUE)
# not really overdispersed
testZeroInflation(simulationOutput_bm, plot = TRUE, alternative = "greater")
## But a bit zeroinflated


# b1 <- lmer(logValue ~  (ESA * ) + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
#              (bio10_4_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePETSD + SnowMonths_cat)^2 + 
#              scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
#              scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
#              ScalePETSD:bio10_12_scaled + ScalePETSD:bio10_15_scaled +
#              (1|file/Study_Name), data = biomass,
#            control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))
# 

biomass_model <- modelSimplificationAIC_glmmTMB(model = b1, itermax = 2e5, evalmax=2e5, dat = biomass)
save(biomass_model, file = file.path(models, "biomassmodel_functionalgroups_correction.rds"))
beep()
# save(biomass_model, file = file.path(models, "biomassmodel_functionalgroups_revised.rds"))
# load(file.path(models, "biomassmodel_functionalgroups.rds"))

simulationOutput_bm <- simulateResiduals(fittedModel = biomass_model, n = 250)
plot(simulationOutput_bm,quantreg = TRUE)
# so ignore the res versus pred plot, but look at others
testDispersion(simulationOutput_bm, alternative = "greater", plot = TRUE)
# not really overdispersed
testZeroInflation(simulationOutput_bm, plot = TRUE, alternative = "greater")
# a bit zero inflated


############### abundance
abundance$logValue <- log(abundance$value + 1)
hist(abundance$logValue)


a1 <- glmmTMB(logValue ~ (ESA * variable) + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
                (bio10_7_scaled + bio10_12_scaled + bio10_15_scaled + SnowMonths_cat + ScalePET)^2 +
                scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
                scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET + 
          (1|file/Study_Name), data = abundance, 
        ziformula = ~ESA + variable + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC +
          bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + bio10_12_scaled + ScalePET,
        control = glmmTMBControl(optCtrl = list(iter.max = 2e5,eval.max=2e5))) #, optimizer ="bobyqa"))
## Checking for 'larger hessians'. Small (~0) numbers give non-postivie Hessian matrix warning
beep()
ff <- fixef(a1)$zi
round(exp(c(intercept=unname(ff[1]),ff[-1]+ff[1])),3)

simulationOutput_am <- simulateResiduals(fittedModel = a1, n = 250)
plot(simulationOutput_am,quantreg = TRUE)
# so ignore the res versus pred plot, but look at others
testDispersion(simulationOutput_am, alternative = "greater", plot = TRUE)
# not really overdispersed
testZeroInflation(simulationOutput_am, plot = TRUE, alternative = "greater")
## Slightly zero inflated, but not really


abundance_model <- modelSimplificationAIC_glmmTMB(model = a1, itermax = 2e5, evalmax=2e5, dat = abundance)
save(abundance_model, file = file.path(models, "abundancemodel_functionalgroups_correction.rds"))
beep()
# save(abundance_model, file = file.path(models, "abundancemodel_functionalgroups_revised.rds"))
# load(file.path(models, "abundancemodel_functionalgroups.rds"))




############## Richness
richness <- m_sites[grep("richness", m_sites$variable),]
richness <- droplevels(richness[which(!(is.na(richness$value))),])
# norichness <- tapply(richness$value, richness$Study_Name, FUN = mean) # often no variation, so use mean
# norichness <- names(norichness[which(is.na(norichness))])
# richness <- droplevels(richness[!(richness$Study_Name %in% norichness),]) # 37980


richness <- droplevels(richness[richness$ESA != "Unknown",]) # 36955


richness <- droplevels(richness[!(is.na(richness$bio10_15)),]) ##  36885
richness <- droplevels(richness[!(is.na(richness$OCFinal)),]) ## 
richness <- droplevels(richness[!(is.na(richness$phFinal)),]) ## 
richness <- droplevels(richness[!(is.na(richness$SnowMonths_cat)),]) ##  
richness <- droplevels(richness[!(is.na(richness$Aridity)),]) ##  36355



table(richness$ESA, richness$variable)

richness_notinclude <- c("Needleleaf deciduous forest", "Tree open", "Sparse vegetation", "Shrub", 
                        "Cropland/Other vegetation mosaic", "Bare area (consolidated", "Wetland/Herbaceous", "Water bodies")
richness <- droplevels(richness[!(richness$ESA %in% richness_notinclude),]) ##   8426 # 36230

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

## This has changed in the correction process

# remove bio 1, 4, 12, petsd, 

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

esaBYfg <- richness %>% # Start by defining the original dataframe, AND THEN...
  group_by(ESA) %>% # Define the grouping variable, AND THEN...
  summarize( # Now you define your summary variables with a name and a function...
    Epi_biomass = sum(value[which(variable == "Epigeic")], na.rm = TRUE),
    Endo_biomass = sum(value[which(variable == "Endogeic")], na.rm = TRUE),
    Ane_biomass = sum(value[which(variable == "Anecic")], na.rm = TRUE),
    EpiEndo_biomass = sum(value[which(variable == "Epi-Endogeic")], na.rm = TRUE),
    Unknown_biomass = sum(value[which(variable == "Unknown")], na.rm = TRUE)
    
  )


as.data.frame(richness %>% group_by( ESA, variable ) %>% summarise( total = sum(value) ))


as.data.frame(esaBYfg)



r1_hurdle <- glmmTMB(value ~  (ESA * variable) + ScaleElevation + 
                       (scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
                       (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                          ScalePET)^2 + 
                       scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                       scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
                       scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
                       (1|file/Study_Name), data = richness,  family = truncated_poisson,# verbose= TRUE,
                     zi = ~ESA + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC +
                       bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                       ScalePET,
                     control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max=2e5)))


richness_test <- richness[richness$variable != "Unknown_richness",]
r1_nounknowns <- glmmTMB(value ~  (ESA * variable) + ScaleElevation + 
                           (scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
                           (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                              ScalePET)^2 + 
                           scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                           scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
                           scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
                           (1|file/Study_Name), data = richness_test,  family = truncated_poisson,# verbose= TRUE,
                         zi = ~ESA + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC +
                           bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                           ScalePET,
                         control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max=2e5)))

## This works....

richness$variable_diff <- as.character(richness$variable)
richness$variable_diff[richness$variable_diff %in% c("Unknown_richness", "EpiEndo_richness")] <- "Other/Unknown"

r1_other <- glmmTMB(value ~  (ESA * variable_diff) + ScaleElevation + 
                           (scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
                           (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                              ScalePET)^2 + 
                           scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                           scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
                           scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
                           (1|file/Study_Name), data = richness,  family = truncated_poisson,# verbose= TRUE,
                         zi = ~ESA + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC +
                           bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                           ScalePET,
                         control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max=2e5)))

## This works....


print("First model done. Now for the simplification process...")
ff <- fixef(r1_hurdle)$zi
round(exp(c(intercept=unname(ff[1]),ff[-1]+ff[1])),3)


# richness_model <- modelSimplificationAIC_glmmTMB(model = r1_hurdle, itermax = 2e5, evalmax=2e5, dat = richness)
richness_model <- modelSimplificationAIC_glmmTMB(model = r1_other, itermax = 2e5, evalmax=2e5, dat = richness)

save(richness_model, file = file.path(data_out, "richnessmodel_functionalgroups_other_correction.rds"))


