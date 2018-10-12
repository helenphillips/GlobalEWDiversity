## To see whether including measured soil properties
## and soil grids makes a difference

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


df_variables_sensitivity <- function(data){
  toMatch <- c("^PHIHOX$",
               "^CLYPPT$", "^SLTPPT$",
               "^CECSOL$", "^ORCDRC$",
               
               "^bio10_1$", "^bio10_4$",
               "^bio10_7$", "^bio10_12$",
               "^bio10_15$",
               
               "^Aridity$", "^PETyr$", "^PET_SD$", "^elevation$")
  
  matches <- unique (grep(paste(toMatch,collapse="|"), 
                          names(data), value=FALSE))
  
  return(matches)
  
}


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


biomass$scaleAridity <- scale(biomass$Aridity)
biomass$ScalePET <- scale(biomass$PETyr)
biomass$ScalePETSD <- scale(biomass$PET_SD)
biomass$ScaleElevation <- scale(biomass$elevation)
## Save the data
# write.csv(biomass, file = file.path(data_out, paste("sitesBiomass_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

## 
ind <- df_variables_sensitivity(biomass)
dat <- biomass[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# Original model:
# bio10_7, bio10_12  bio10_15, CECSOL, 
# elevation, PETyr,phFinal,ClayFinal,SiltFinal,OCFinal

# With only soilgrids data:
# bio10_12 bio10_15  PHIHOX    CLYPPT   SLTPPT    CECSOL
# ORCDRC   elevation PET_SD   


b1 <- lmer(logBiomass ~  ESA + ScaleElevation + 
             (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
             (bio10_12_scaled  + bio10_15_scaled + ScalePETSD + SnowMonths_cat)^2 + 
             scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
             ScalePETSD:bio10_12_scaled + ScalePETSD:bio10_15_scaled +
             (1|file/Study_Name), data = biomass,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))


biomass_model_SG <- modelSimplificationAIC(model = b1, data = biomass, optimizer = "bobyqa", Iters = 2e5)
save(biomass_model_SG, file = file.path(models, "biomassmodel_SoilGrids.rds"))

# Load the original model
load(file.path(models, "biomassmodel_full.rds"))

summary(biomass_model_SG)


new <- data.frame(fixef(biomass_model_SG))
new$var <- rownames(new)

old <- data.frame(fixef(biomass_model))
old$var <- rownames(old)

biomass_all <- merge(new, old, by = "var", all = TRUE)
biomass_all[is.na(biomass_all)] <- 10

plot(biomass_all$fixef.biomass_model_SG. ~ biomass_all$fixef.biomass_model., 
     pch = 19, ylab = "Only SoilGrids data", xlab = "Mixture of data")
abline(0,1)


#################################################
# 6. Abundance
#################################################
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

abundance$scaleAridity <- scale(abundance$Aridity)
abundance$ScalePET <- scale(abundance$PETyr)
abundance$ScalePETSD <- scale(abundance$PET_SD)
abundance$ScaleElevation  <- scale(abundance$elevation)
## Save the data
# write.csv(abundance, file = file.path(data_out, paste("sitesAbundance_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

ind <- df_variables_sensitivity(abundance)
dat <- abundance[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# Original model
# bio10_7, bio10_15,CECSOL,elevation,Aridity,PETyr,  phFinal,ClayFinal,SiltFinal, OCFinal   

## With only soil grids
# bio10_7   bio10_15  PHIHOX    CLYPPT  SLTPPT   CECSOL   ORCDRC   elevation Aridity   PETyr 

a1 <- lmer(logAbundance ~  ESA + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
             (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                ScalePET)^2 +
             scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
             scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET + 
             scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity + 
             (1|file/Study_Name), data = abundance,
           control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

abundance_model_SG <- modelSimplificationAIC(model = a1, data = abundance, optimizer = "bobyqa", Iters = 2e5)
save(abundance_model_SG, file = file.path(models, "abundancemodel_SoilGrids.rds"))


load(file.path(models, "abundancemodel_full.rds"))



new <- data.frame(fixef(abundance_model_SG))
new$var <- rownames(new)

old <- data.frame(fixef(abundance_model))
old$var <- rownames(old)

abundance_all <- merge(new, old, by = "var", all = TRUE)
abundance_all[is.na(abundance_all)] <- 10

plot(abundance_all$fixef.abundance_model_SG. ~ abundance_all$fixef.abundance_model., 
     pch = 19, ylab = "Only SoilGrids data", xlab = "Mixture of data")
abline(0,1)
