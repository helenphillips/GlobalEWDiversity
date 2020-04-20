args <- commandArgs(trailingOnly = TRUE)

data_out <- args[1] # output directory
data_in <- args[2] # data folder
date <- args[3] # date for dataset
functions <- args[4] # folder containing my functions
##############
print(functions)
##### Functions
library(glmmTMB)
source(file.path(functions, "ModelSimplification.R"))
source(file.path(functions, "FormatData.R"))

##### Data
## Script 4 does the analysis to find out which variables are needed in the model 
# Variables have also been scaled
print(paste("Loading the file: ", file.path(data_in, paste("sitesRichness_", date, ".csv", sep = ""))))
richness <- read.csv(file = file.path(data_in, paste("sitesRichness_", date, ".csv", sep = "")))

##### Modelling
print("Doing the first model")

r1_hurdle <- glmmTMB(SpeciesRichness ~  ESA + scaleElevation + 
                       (scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
                       (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                          ScalePET)^2 + 
                       scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                       scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET +
                       scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity +
                       (1|file/Study_Name), data = richness,  family = truncated_poisson,# verbose= TRUE,
                     zi = ~ESA + scaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC +
                       bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                       ScalePET,
                     control = glmmTMBControl(optCtrl = list(iter.max = 2e5, eval.max=2e5)))

print("First model done. Now for the simplification process...")


richness_model <- modelSimplificationAIC_glmmTMB(model = r1_hurdle, itermax = 2e5, evalmax=2e5, dat = richness)

save(richness_model, file = file.path(data_out, "richnessmodel_correction.rds"))

print("Done!")