
if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 0. Loading libraries
#################################################
library(maptools)
library(maps)
library(lme4)
library(MCMCglmm)
library(car)
library(DHARMa)
library(reshape)
source("Functions/FormatData.R")
source("Functions/CorvifVariablePicker.R")
source("Functions/lme4_ModellingFunctions.R")
source("Functions/ModelSimplification.R")
source("MEE3_1_sm_Appendix_S1/HighstatLib.R")


#################################################
# 1. Create folders
#################################################
models <- "Models"
if(!dir.exists("15_Data")){
  dir.create("15_Data")
}

data_out <- "15_Data"

#################################################
# 2. Loading in species data
#################################################

data_in <-"10_Data"

files <- list.files(file.path(data_in))
files <- files[grep("SiteswithFunctionalGroups_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

rm(files)
rm(date)

#################################################
# 5. Load in data
#################################################

sites <- read.csv(file.path(data_in, loadin))
# sites <- read.csv("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\3_Data\\Sites_2017-11-09.csv")
rm(loadin)

sites <- SiteLevels(sites) ## relevels all land use/habitat variables


table(sites$FGRichness)
##############################################################
##
#############################################################

## 
sites <- droplevels(sites[!(is.na(sites$bio10_15)),]) ## 
sites <- droplevels(sites[!(is.na(sites$OCFinal)),]) ## 
sites <- droplevels(sites[!(is.na(sites$phFinal)),]) ## 
sites <- droplevels(sites[!(is.na(sites$scaleAridity)),]) ## 
sites <- droplevels(sites[!(is.na(sites$SnowMonths_cat)),]) ## 

sites <- droplevels(sites[sites$ESA != "Unknown",]) # 


table(sites$ESA)
sites_notinclude <- c("Needleleaf deciduous forest", "Tree open",
                         "Sparse vegetation",  "Urban",
                         "Bare area (consolidated", "Paddy field", "Wetland/Herbaceous", "Water bodies")

sites <- droplevels(sites[!(sites$ESA %in% sites_notinclude),]) ##   5363
summary(sites$phFinal)
sites$scalePH <- scale(sites$phFinal)
sites$scaleCLYPPT <- scale(sites$ClayFinal)
sites$scaleSLTPPT <- scale(sites$SiltFinal)
sites$scaleCECSOL <- scale(sites$CECSOL)
sites$scaleORCDRC <- scale(sites$OCFinal)

sites$bio10_1_scaled <- scale(sites$bio10_1)
sites$bio10_4_scaled <- scale(sites$bio10_4)
sites$bio10_7_scaled <- scale(sites$bio10_7)
sites$bio10_12_scaled <- scale(sites$bio10_12)
sites$bio10_15_scaled <- scale(sites$bio10_15)

sites$scaleAridity <- scale(sites$Aridity)
sites$ScalePET <- scale(sites$PETyr)
sites$ScalePETSD <- scale(sites$PET_SD)
sites$ScaleElevation <- scale(sites$elevation)

## Save the data
write.csv(sites, file = file.path(data_out, paste("sites+FGRichness_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

ind <- df_variables(sites)
dat <- sites[,c(ind)]
cor <- findVariables(dat, VIFThreshold = 3)

# "bio10_7"   "bio10_15"  "CECSOL"  "elevation"  "Aridity"   "PETyr"     "phFinal"  
# "ClayFinal" "SiltFinal" "OCFinal"  


# Run the model of the cluster


fg1 <- glmer(FGRichness ~  ESA + ScaleElevation + 
               (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
              (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                 ScalePET)^2 + 
              scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
              scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
              scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
              (1|file/Study_Name), data = sites, family = poisson,
            control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))


summary(fg1)
save(fg1, file = file.path(models, "functionalrichnessmodel_initialmodel.rds"))
load(file.path(models, "functionalrichnessmodel_initialmodel.rds"))

#####
# Looking at coefficients of final model
load(file.path(models, "fgrichnessmodel.rds"))


summary(fgrichness_model)

## Can use it to predict the values we already
## This will help me see if the estimates are reasonable

source("Functions/Plots.R")
sites <- droplevels(sites[sites$ESA != "Bare area (unconsolidated",])
sites$Predicted <- (predict(fgrichness_model, sites, re.form = NA))

summary(exp(predictedSites[,136]))
plot(exp(sites$Predicted) ~ jitter(sites$FGRichness))
abline(0, 1)
## Still underpredicts at the higher FGrichness....