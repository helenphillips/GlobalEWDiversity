if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}
#################################################
# 1. Libraries
#################################################
library(lme4)

createSplits <- function(dat, kfold = 10){
  rows <- nrow(dat)
  rows <- sample(rows, size = length(1:rows))
  
  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  
  splits <- chunk(rows, kfold)
  
  return(splits)
}
#################################################
# 2. Load data
#################################################



#################################################
# 3. Create directories
#################################################

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"


#################################################
# 4. Load in models
#################################################

models <- "Models"


load(file.path(models, "richnessmodel_full.rds"))
load(file.path(models, "biomassmodel_full.rds"))
load(file.path(models, "abundancemodel_full.rds"))

#################################################
# 5. Species Richness
################################################

spRData <- richness_model@frame

spRData$Predicted <- exp(predict(richness_model, spRData, re.form = NA))

opaqueBlack <- rgb(red = 0, green = 0, blue = 0, alpha = 0.3,maxColorValue = 1)
plot(spRData$Predicted ~ jitter(spRData$SpeciesRichness), ylim = c(0, 15), pch = 19, col = opaqueBlack)
abline(0, 1) 


#################################################
# 5. Abundance
################################################

abundanceData <- abundance_model@frame

abundanceData$Predicted <- (predict(abundance_model, abundanceData, re.form = NA))

opaqueBlack <- rgb(red = 0, green = 0, blue = 0, alpha = 0.3,maxColorValue = 1)
plot(abundanceData$Predicted ~ jitter(abundanceData$logAbundance), ylim = c(0, 8), pch = 19, col = opaqueBlack)
abline(0, 1) 

########
# K-Fold Cross validation
########
  k_fold <- 10
splits <- createSplits(abundanceData, kfold = k_fold)

predictedData <- list()
for(k in 1:k_fold){
  
  rows <- as.vector(unlist(splits[k]))
  testData <- abundanceData[rows,]
  bankData <- abundanceData[-rows,]
  
  mod <- lmer(logAbundance ~ scalePH + scaleCLYPPT + scaleSLTPPT + scaleCECSOL +  
                scaleORCDRC + bio10_1_scaled + bio10_15_scaled + SnowMonths_cat +  
                scaleAridity + ScalePETSD + scalePH:scaleCLYPPT + scalePH:scaleCECSOL +  
                scaleCLYPPT:scaleCECSOL + scaleSLTPPT:scaleCECSOL + scaleCECSOL:scaleORCDRC +  
                bio10_1_scaled:bio10_15_scaled + bio10_1_scaled:SnowMonths_cat +  
                bio10_1_scaled:scaleAridity + bio10_1_scaled:ScalePETSD +  
                bio10_15_scaled:SnowMonths_cat + SnowMonths_cat:ScalePETSD +  
                scaleCLYPPT:bio10_15_scaled + scaleCLYPPT:ScalePETSD + ESA +  
                (1 | file/Study_Name), data = bankData, 
              control = lmerControl(optCtrl = list(maxfun = 2e5), optimizer ="bobyqa"))
  
  testData$Predicted <- (predict(mod, testData, re.form = NA))
  
  predictedData[[k]] <- data.frame(observed = testData$logAbundance, predicted = testData$Predicted)
  
}

df <- do.call("rbind", predictedData)
plot(df$predicted ~ df$observed)
abline(0, 1) 
