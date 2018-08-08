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

epi_richness <- fg_richness[grep("Epi", fg_richness$variable),]
epi_richness <- scaleVariables(epi_richness)

##### Modelling
# Scaling variables in this script, but using VIF results from previous
print("Doing the first model")

epi_r1 <- glmer(value ~  ESA + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleORCDRC + scaleCECSOL)^2 +
                  (bio10_12_scaled + bio10_15_scaled + SnowMonths_cat + ScalePET + ScalePETSD)^2 + 
                  scaleCLYPPT:bio10_12_scaled + scaleSLTPPT:bio10_12_scaled +
                  scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                  scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
                  scaleCLYPPT:ScalePETSD + scaleSLTPPT:ScalePETSD +
                  
                  (1|file/Study_Name), data = epi_richness, family = poisson,
                control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

print("First model done. Now for the simplification process...")


epi_richness_model <- modelSimplificationAIC(model = epi_r1, data = epi_richness, optimizer = "bobyqa", Iters = 2e5)
save(epi_richness_model, file = file.path(data_out, "richnessmodel_epifunctionalgroups.rds"))

print("Done!")