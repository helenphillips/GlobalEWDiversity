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
## Script 4 does the analysis to find out which variables are needed in the model 
# Variables have also been scaled
print(paste("Loading the file: ", file.path(data_in, paste("sitesRichness_", date, ".csv", sep = ""))))
richness <- read.csv(file = file.path(data_in, paste("sitesRichness_", date, ".csv", sep = "")))

##### Modelling
print("Doing the first model")

r1 <- glmer(SpeciesRichness ~  ESA + scaleElevation + 
              (scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
              (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                 ScalePET)^2 + 
              scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
              scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
              scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
              (1|file/Study_Name), data = richness, family = poisson,
            control = glmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))

print("First model done. Now for the simplification process...")


richness_model <- modelSimplificationAIC(model = r1, data = richness, optimizer = "bobyqa", Iters = 2e5)
save(richness_model, file = file.path(data_out, "richnessmodel_correction.rds"))

print("Done!")