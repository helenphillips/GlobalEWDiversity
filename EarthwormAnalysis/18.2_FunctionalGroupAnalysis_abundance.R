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
print(paste("Loading the file: ", file.path(data_in, paste("sitesFGAbundance_", date, ".csv", sep = ""))))
abundance <- read.csv(file = file.path(data_in, paste("sitesFGAbundance_", date, ".csv", sep = "")))

##### Modelling

print("Doing the first model")

a1 <- glmmTMB(logValue ~ (ESA * variable) + ScaleElevation + (scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC)^2 +
                (bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + 
                   ScalePET)^2 +
                scaleCLYPPT:bio10_15_scaled + scaleSLTPPT:bio10_15_scaled +
                scaleCLYPPT:ScalePET + scaleSLTPPT:ScalePET + 
                scaleCLYPPT:scaleAridity + scaleSLTPPT:scaleAridity + 
          (1|file/Study_Name), data = abundance, 
        ziformula = ~ESA + variable + ScaleElevation + scalePH  + scaleCLYPPT + scaleSLTPPT + scaleCECSOL + scaleORCDRC +
          bio10_7_scaled + bio10_15_scaled + SnowMonths_cat + scaleAridity + ScalePET,
        control = glmmTMBControl(optCtrl = list(iter.max = 2e5,eval.max=2e5))) #, optimizer ="bobyqa"))


print("First model done. Now for the simplification process...")

abundance_model <- modelSimplificationAIC_glmmTMB(model = a1, itermax = 2e5, evalmax=2e5, dat = abundance)
save(abundance_model, file = file.path(data_out, "abundancemodel_functionalgroups.rds"))

print("Done!")
