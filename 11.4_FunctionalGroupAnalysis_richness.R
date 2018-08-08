args <- commandArgs(trailingOnly = TRUE)

data_out <- args[1] # output directory
data_in <- args[2] # data folder
date <- args[3] # date for dataset
functions <- args[4] # folder containing my functions
##############

##### Functions

source(file.path(functions, "ModelSimplification.R"))
source(file.path(functions, "FormatData.R"))

##### Data
## Script 11 does the analysis to find out which variables are needed in the model 
# Variables have also been scaled
print(paste("Loading the file: ", file.path(data_in, paste("sitesFGRichness_", date, ".csv", sep = ""))))
fg_richness <- read.csv(file = file.path(data_in, paste("sitesFGRichness_", date, ".csv", sep = "")))

ane_richness <- fg_richness[grep("Ane", fg_richness$variable),]
ane_richness <- scaleVariables(ane_richness)

##### Modelling
# Scaling variables in this script, but using VIF results from previous
print("Doing the first model")

ane_r1 <- glmer(value ~  ESA + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                  (bio10_1_scaled + bio10_12_scaled + bio10_15_scaled + SnowMonths_cat + ScalePETSD)^2 + 
                  scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                  scaleCLYPPT:ScalePETSD + scaleSLTPPT:ScalePETSD +
                  scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
                  
                  (1|file/Study_Name), data = ane_richness, family = poisson,
                control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

print("First model done. Now for the simplification process...")


ane_richness_model <- modelSimplificationAIC(model = ane_r1, data = ane_richness, optimizer = "bobyqa", Iters = 2e5)
save(ane_richness_model, file = file.path(data_out, "richnessmodel_anefunctionalgroups.rds"))

print("Done!")