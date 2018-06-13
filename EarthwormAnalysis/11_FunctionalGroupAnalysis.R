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
# 6. Split into (currently) two diversity measures
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
write.csv(biomass, file = file.path(data_out, paste("sitesFGBiomass_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)


corvif(data.frame(biomass$bio10_1,biomass$bio10_4,biomass$bio10_7,biomass$bio10_12,biomass$bio10_15, 
                  biomass$Aridity, biomass$PETyr, biomass$PET_SD,
                  biomass$phFinal, biomass$ClayFinal, biomass$SiltFinal, biomass$OCFinal, biomass$CECSOL))
# Remove bio 7
corvif(data.frame(biomass$bio10_1,biomass$bio10_4,biomass$bio10_12,biomass$bio10_15, 
                  biomass$Aridity, biomass$PETyr, biomass$PET_SD,
                  biomass$phFinal, biomass$ClayFinal, biomass$SiltFinal, biomass$OCFinal, biomass$CECSOL))
# REmove 1
corvif(data.frame(biomass$bio10_4,biomass$bio10_12,biomass$bio10_15, 
                  biomass$Aridity, biomass$PETyr, biomass$PET_SD,
                  biomass$phFinal, biomass$ClayFinal, biomass$SiltFinal, biomass$OCFinal, biomass$CECSOL))
# Remove petyr
corvif(data.frame(biomass$bio10_4,biomass$bio10_12,biomass$bio10_15, 
                  biomass$Aridity,biomass$PET_SD,
                  biomass$phFinal, biomass$ClayFinal, biomass$SiltFinal, biomass$OCFinal, biomass$CECSOL))
# Remove Aridity
corvif(data.frame(biomass$bio10_4,biomass$bio10_12,biomass$bio10_15, 
                  biomass$PET_SD,
                  biomass$phFinal, biomass$ClayFinal, biomass$SiltFinal, biomass$OCFinal, biomass$CECSOL))
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

write.csv(abundance, file = file.path(data_out, paste("sitesFGAbundance_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

corvif(data.frame(abundance$bio10_1,abundance$bio10_4,abundance$bio10_7,abundance$bio10_12,abundance$bio10_15, 
                  abundance$Aridity, abundance$PETyr, abundance$PET_SD,
                  abundance$phFinal, abundance$ClayFinal, abundance$SiltFinal, abundance$OCFinal, abundance$CECSOL))
# Remove bio 7
corvif(data.frame(abundance$bio10_1,abundance$bio10_4,abundance$bio10_12,abundance$bio10_15, 
                  abundance$Aridity, abundance$PETyr, abundance$PET_SD,
                  abundance$phFinal, abundance$ClayFinal, abundance$SiltFinal, abundance$OCFinal, abundance$CECSOL))
# Remove 1
corvif(data.frame(abundance$bio10_4,abundance$bio10_12,abundance$bio10_15, 
                  abundance$Aridity, abundance$PETyr, abundance$PET_SD,
                  abundance$phFinal, abundance$ClayFinal, abundance$SiltFinal, abundance$OCFinal, abundance$CECSOL))
# Remove 12
corvif(data.frame(abundance$bio10_4,abundance$bio10_15, 
                  abundance$Aridity, abundance$PETyr, abundance$PET_SD,
                  abundance$phFinal, abundance$ClayFinal, abundance$SiltFinal, abundance$OCFinal, abundance$CECSOL))
# Remove pet_sd
corvif(data.frame(abundance$bio10_4,abundance$bio10_15, 
                  abundance$Aridity, abundance$PETyr,
                  abundance$phFinal, abundance$ClayFinal, abundance$SiltFinal, abundance$OCFinal, abundance$CECSOL))
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
