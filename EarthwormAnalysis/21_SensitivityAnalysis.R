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

data_in <-"3.5_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

rm(files)
rm(date)

sites <- read.csv(file.path(data_in, loadin))

data_in <-"4_Data"
files <- list.files(file.path(data_in))
files <- files[grep("sitesBiomass_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

biomass <- read.csv(file.path(data_in, loadin))


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

## Which have "bad" predictions

badPreds <- droplevels(spRData[(spRData$SpeciesRichness) > 6,])
levels(badPreds$file)

dat <- droplevels(sites[sites$file %in% badPreds$file,])
levels(dat$Country)

summary(spRData)
summary(badPreds)

library(maps)
library(maptools)
coord<-aggregate(cbind(dat$Longitude__Decimal_Degrees, dat$Latitude__decimal_degrees), list(dat$Study_Name), mean)

coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")
dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")
mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)


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


### Biomass
dat <- biomass_model@frame
dat$resids <- residuals(biomass_model)
dat$predicted <- (predict(biomass_model, dat, re.form = NA))
plot(dat$predicted ~ dat$logBiomass)
abline(0, 1) 


