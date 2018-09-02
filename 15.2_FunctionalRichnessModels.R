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
print(paste("Loading the file: ", file.path(data_in, paste("sites+FGRichness_", date, ".csv", sep = ""))))
fg_richness <- read.csv(file = file.path(data_in, paste("sites+FGRichness_", date, ".csv", sep = "")))

##### Modelling
# Scaling variables in this script, but using VIF results from previous
print("Doing the first model")

fg1 <- glmer(FGRichness ~  ESA + (scalePH  + 
                                    scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
               (bio10_1_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                  ScalePETSD)^2 + 
               scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
               scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
               scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
               (1|file/Study_Name), data = fg_richness, family = poisson,
             control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

print("First model done. Now for the simplification process...")


fgrichness_model <- modelSimplificationAIC(model = fg1, data = fg_richness, optimizer = "bobyqa", Iters = 2e5)
save(fgrichness_model, file = file.path(data_out, "fgrichnessmodel.rds"))

print("Done!")