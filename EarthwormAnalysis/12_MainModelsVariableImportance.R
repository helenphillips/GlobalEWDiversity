########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################

library(lme4)
library(randomForest)
library(ggplot2)
library(reshape)
library(viridis)

source("Functions/FormatData.R")
source("Functions/lme4_ModellingFunctions.R")
source("Functions/ModelSimplification.R")
source("MEE3_1_sm_Appendix_S1/HighstatLib.R")
# source("Functions/HypothesisTesting.R")
source("Functions/RF_VariableSetImportance.R")
source("Functions/Plots.R")
#################################################
# 2. Variables
#################################################


data_in <- "8_Data"
models <- "Models"
figures <- "Figures"


circleSize <- function(dat){
  dat$size <- 9
  onepercent <- min(dat$delta, na.rm = TRUE) / 100
  dat$percentChange <- dat$delta / onepercent
  dat$size <- dat$size * ((100 - dat$percentChange) / 100)
  dat$size <- dat$size + 1
  return(dat)
}
#################################################
# 3. Load data and models
#################################################
## Richness 
files <- list.files(file.path(data_in))
files <- files[grep("sitesRichness_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
richness <- read.csv(file.path(data_in, loadin))

load(file.path(models, "richnessmodel_revised.rds"))

## Biomass
files <- list.files(file.path(data_in))
files <- files[grep("sitesBiomass_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
biomass <- read.csv(file.path(data_in, loadin))

load(file.path(models, "biomassmodel_full_revised.rds"))

## Abundance
files <- list.files(file.path(data_in))
files <- files[grep("sitesAbundance_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
abundance <- read.csv(file.path(data_in, loadin))

load(file.path(models, "abundancemodel_full_revised.rds"))




#################################################
# With Different Groups
##################################################

# ABUNDANCE


AbundmainEffects <- c("scalePH" , "scaleCLYPPT" ,"scaleSLTPPT" , "scaleCECSOL" , 
                      "scaleORCDRC" , "bio10_7_scaled" , "bio10_15_scaled" , "SnowMonths_cat" ,
                      "scaleAridity" , "ScalePET", "ESA", "ScaleElevation")
# 12
abundance_rf <- randomForest(y = abundance$logAbundance, x = abundance[,names(abundance) %in% AbundmainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)


# my_varimp <- importance(abundance_rf, scale=FALSE)

ESA <- "ESA"
Elevation <- "ScaleElevation"
Temperature <- c("bio10_7_scaled", "ScalePET")
Precip <- c("bio10_15_scaled" , "SnowMonths_cat", "scaleAridity")
Soil <- c("scalePH", "scaleCLYPPT", "scaleCECSOL", 
          "scaleSLTPPT", "scaleORCDRC")
WaterRetention <- c("scaleCLYPPT", "bio10_15_scaled", "ScalePET", "scaleAridity")

groups <- list(
  ESA = ESA,
  Elevation = Elevation, 
  Temperature = Temperature,
  Precip = Precip,
  Soil = Soil,
  WaterRetention = WaterRetention
)
abundance_import_split <- group.importance(abundance_rf, groups) 


# RICHNESS
richness_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL",  
                          "scaleORCDRC", "bio10_7_scaled", "bio10_15_scaled", "SnowMonths_cat",  
                          "scaleAridity", "ScalePET", "ESA", "scaleElevation")
# 12
spR_rf <- randomForest(y = richness$SpeciesRichness, x = richness[,names(richness) %in% richness_mainEffects], 
                       ntree=501, importance=TRUE, proximity = TRUE)



ESA <- "ESA"
Elevation <- "scaleElevation"
Temperature <- c("bio10_7_scaled", "ScalePET")
Precip <- c("bio10_15_scaled", "SnowMonths_cat", "scaleAridity")
Soil <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL", "scaleORCDRC")
WaterRetention <- c("scaleCLYPPT","scaleSLTPPT", "bio10_15_scaled", "scaleAridity", "ScalePET")

groups <- list(
  ESA = ESA,
  Elevation = Elevation,
  Temperature = Temperature,
  Precip = Precip, 
  Soil = Soil,
  WaterRetention = WaterRetention
)
richness_import_split <- group.importance(spR_rf, groups) 

# BIOMASS
biomass_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "scaleCECSOL", 
                         "bio10_7_scaled",  "bio10_12_scaled", 
                         "bio10_15_scaled", "ScalePET", "SnowMonths_cat", "ESA")
# 11
bioM_rf <- randomForest(y = biomass$logBiomass, x = biomass[,names(biomass) %in% biomass_mainEffects], 
                        ntree=501, importance=TRUE, proximity = TRUE)

ESA <- "ESA"
Temperature <- c("bio10_7_scaled", "ScalePET")
Precip <- c("bio10_12_scaled", "bio10_15_scaled", "SnowMonths_cat")
Soil <- c("scalePH","scaleORCDRC", "scaleCECSOL","scaleCLYPPT", "scaleSLTPPT")
WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_12_scaled", "bio10_15_scaled")

groups <- list(
  ESA = ESA,
  Temperature = Temperature,
  Precip = Precip,
  Soil = Soil,
  WaterRetention = WaterRetention
)
biomass_import_split <- group.importance(bioM_rf, groups) 


######################
# Ordering
#######################
abundance_import_split
abundance_order <- c(1, 3, 4, 6, 2, 5)
richness_import_split
richness_order <- c(4, 5, 6, 1, 2, 3)
biomass_import_split
biomass_order <- c(5,NA, 6,4,3, 2)

abundance_import_delta <- c(abundance_import_split - max(abundance_import_split, na.rm = TRUE))
richness_import_delta <- c(richness_import_split - max(richness_import_split, na.rm = TRUE))
biomass_import_delta <- c(biomass_import_split - max(biomass_import_split, na.rm = TRUE))
biomass_import_delta <- c(biomass_import_delta[1], NA, biomass_import_delta[2:5])


d <- a <- matrix(rep(NA, length = 6*3), nrow = 3, ncol = 6)
colnames(d) <- colnames(a) <- c("Habitat Cover","Elevation", "Temperature", "Precipitation", "Soil", "Water Retention")
a[1,] <- richness_order
a[2,] <- abundance_order
a[3,] <- biomass_order

## 4 is most important, 1 is least important
rownames(d) <- rownames(a) <- c("Species Richness", "Abundance", "Biomass")


d[1,] <- richness_import_delta
d[2,] <- abundance_import_delta
d[3,] <- biomass_import_delta


d <- round(d, digits = 0)

dat <- melt(a)


dat$X1 <- factor(dat$X1, levels = c( "Biomass", "Abundance","Species Richness"))
dat$X2 <- factor(dat$X2, levels = c("Habitat Cover","Elevation", "Soil","Precipitation", "Temperature","Water Retention"))



delta <- melt(d)

alldat <- cbind(dat, delta$value)
names(alldat)[names(alldat) == "delta$value"] <- "delta"

ord <- c(1, 2, 4, 5, 3, 6)

spR <- alldat[alldat$X1 == "Species Richness",]
spR$y <- 3
spR$x <- ord
abund <- alldat[alldat$X1 == "Abundance",]
abund$y <- 2
abund$x <- ord
bmass<- alldat[alldat$X1 == "Biomass",]
bmass$y <- 1
bmass$x <- ord



spR <- circleSize(spR)
bmass<- circleSize(bmass)
abund<- circleSize(abund)


all_dat <- rbind(spR, abund, bmass)

labs <- c("Habitat Cover","Elevation","Soil","Precipitation","Temperature","Water\nRetention")

jpeg(file = file.path(figures, "variableImportance_splitGroups_circles.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
par(mar = c(3, 9.5, 1, 5))
plot(-1e+05, -1e+05, ylim = c(0, 4), xlim = c(0.5, 6.5),  
     ylab = "", xlab = "",  xaxt='n', axes = FALSE)
axis(side = 2, cex.axis = 1, labels = levels(all_dat$X1), 
     at = c(1:3), las = 2)
points(all_dat$x, all_dat$y, pch = 19, cex = all_dat$size, ylim = c(0, 5))
axis(side=1, at = 1:6, labels = labs, las=1, cex.axis = 1, padj=1, mgp = c(3, 0, 0)) 
# padj put all labels on the same line, and mgp puts labels closer to the axis
mtext("Model", side = 2, line = 6, cex = 2)
dev.off()


########################################################################
##
########################################################################

# Better to use MSE than NodeInpurity
# https://stats.stackexchange.com/questions/162465/in-a-random-forest-is-larger-incmse-better-or-worse


rich <- importance(spR_rf, scale = FALSE)[,1]
abund <- importance(abundance_rf)[,1]
bio <- importance(bioM_rf)[,1]


allVars <- unique(c(names(rich), names(abund),names(bio)))
allVars[which(allVars == "ScaleElevation")] <- "scaleElevation"
allVars <- unique(allVars)
allVars <- as.factor(allVars)
allVars <- droplevels(allVars)
allVars <- factor(allVars, levels = c(
  "ESA",
  "scaleElevation",
  #soil
  "scaleCECSOL",
  "scaleCLYPPT",
  "scaleORCDRC","scalePH",
  "scaleSLTPPT",
  #temp
  "bio10_7_scaled","ScalePET",
  # precip
  "bio10_15_scaled","SnowMonths_cat",  "scaleAridity","bio10_12_scaled"
))



getPlottingDf <- function(vect, allVars){
  
  df <- data.frame(delta = (vect - vect[which(vect == max(vect))]))
  df <- circleSize(df)
  
  
  for(r in 1:nrow(df)){
    if(rownames(df)[r] == "ScaleElevation"){
      df$x[r] <- 2
    }else {
      df$x[r] <- grep(rownames(df)[r], levels(allVars))
    }
  }
  
  return(df)
}

rich <- getPlottingDf(rich, allVars)
rich$y <- 3
abund <- getPlottingDf(abund, allVars)
abund$y <- 2
bio <- getPlottingDf(bio, allVars)
bio$y <- 1


allDat <- rbind(rich, abund, bio)

labs <- c( "Habitat Cover","Elevation",  "CEC"    ,
          "Clay","Org C","PH"        ,
          "Silt","bio10_7" , "PET"       ,
          "bio10_15", "Snow Months","Aridity Index"   ,
          "bio10_12")

jpeg(file = file.path(figures, "variableImportanceMSE_splitGroups_circles.jpg"), quality = 100, res = 200, width = 2000, height = 1000)

par(mar = c(10, 7, 1, 1))
plot(-1e+05, -1e+05, ylim = c(0, 4), xlim = c(0, length(allVars)),
     ylab = "", xlab = "",  xaxt='n', axes = FALSE)
axis(side = 2, cex.axis = 1, labels = c( "Biomass", "Abundance", "Species Richness"), 
     at = c(1:3), las = 2)
points(allDat$x, allDat$y, pch = 19, cex = allDat$size, ylim = c(0, 5))
axis(side=1, at = 1:length(allVars), labels = labs, las=2, cex.axis = 1) 

axis(1,at=c(2.8, 3, 4, 5, 6, 7, 7.2),col="blue",line=8,tick=T,labels=rep("",7),lwd=2,lwd.ticks=0)
axis(1,at=c(7.8, 8, 9, 9.2),col="blue",line=8,tick=T,labels=rep("",4),lwd=2,lwd.ticks=0)
axis(1,at=c(9.8, 10, 11, 12, 13, 13.2),col="blue",line=8,tick=T,labels=rep("",6),lwd=2,lwd.ticks=0)
mtext("Soil", side = 1, line = 7, at = 5, cex = 1.2)
mtext("Temperature", side = 1, line = 7, at = 8.5, cex = 1.2)
mtext("Precipitation", side = 1, line = 7, at = 11.5, cex = 1.2)

mtext("Model", side = 2, line = 6, cex = 2)
dev.off()
