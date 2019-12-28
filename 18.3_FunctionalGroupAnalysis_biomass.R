args <- commandArgs(trailingOnly = TRUE)

data_out <- args[1] # output directory
data_in <- args[2] # data folder
date <- args[3] # date for dataset
functions <- args[4] # folder containing my functions
##############

##### Functions

source(file.path(functions, "ModelSimplification.R"))

##### Data
## Script 11 does the analysis to find out which variables are needed in the model 
# Variables have also been scaled
print(paste("Loading the file: ", file.path(data_in, paste("sitesFGBiomass_", date, ".csv", sep = ""))))
biomass <- read.csv(file = file.path(data_in, paste("sitesFGBiomass_", date, ".csv", sep = "")))

##### Modelling

print("Doing the first model")

biomass$logValue <- log(biomass$value + 1)
hist(biomass$logValue)

b1 <- glmmTMB(logValue ~  (ESA * variable) + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                (bio10_7_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePET + SnowMonths_cat)^2 + 
                scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
                scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                (1|file/Study_Name), data = biomass, 
              ziformula = ~ESA + variable + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL +
                bio10_7_scaled + bio10_12_scaled  + bio10_15_scaled + ScalePET + SnowMonths_cat,
              control = glmmTMBControl(optCtrl = list(iter.max = 2e5,eval.max=2e5))) #, optimizer ="bobyqa"))


print("First model done. Now for the simplification process...")

biomass_model <- modelSimplificationAIC_glmmTMB(model = b1, itermax = 2e5, evalmax=2e5, dat = biomass)


save(biomass_model, file = file.path(data_out, "biomassmodel_functionalgroups.rds"))

print("Done!")
