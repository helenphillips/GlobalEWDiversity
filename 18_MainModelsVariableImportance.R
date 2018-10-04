########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
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


data_in <- "4_Data"
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

load(file.path(models, "richnessmodel.rds"))

## Biomass
files <- list.files(file.path(data_in))
files <- files[grep("sitesBiomass_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
biomass <- read.csv(file.path(data_in, loadin))

load(file.path(models, "biomassmodel_full.rds"))

## Abundance
files <- list.files(file.path(data_in))
files <- files[grep("sitesAbundance_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
abundance <- read.csv(file.path(data_in, loadin))

load(file.path(models, "abundancemodel_full.rds"))



## Functional Richness
data_in <- "15_Data"

files <- list.files(file.path(data_in))
files <- files[grep("sites\\+FGRichness_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
fgrichness <- read.csv(file.path(data_in, loadin))

load(file.path(models, "fgrichnessmodel.rds"))





#################################################
# With Different Groups
##################################################

# ABUNDANCE


AbundmainEffects <- c("scalePH" , "scaleCLYPPT" ,"scaleSLTPPT" , "scaleCECSOL" , 
                      "scaleORCDRC" , "bio10_7_scaled" , "bio10_15_scaled" , "SnowMonths_cat" ,
                      "scaleAridity" , "ScalePET", "ESA", "ScaleElevation")

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


## Functional Richness

# 
fgrichness_mainEffects <- c("scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "bio10_7_scaled", "bio10_15_scaled",  
                              "SnowMonths_cat", "scaleAridity", "ScalePET", "ScaleElevation", "scalePH")

fgrichness <- droplevels(fgrichness[!(is.na(fgrichness$FGRichness)),])

fgR_rf <- randomForest(y = fgrichness$FGRichness, x = fgrichness[,names(fgrichness) %in% fgrichness_mainEffects], 
                       ntree=501, importance=TRUE, proximity = TRUE)



# ESA <- "ESA"
Elevation <- "ScaleElevation"
Temperature <- c("bio10_7_scaled", "ScalePET")
Precip <- c("bio10_15_scaled", "SnowMonths_cat", "scaleAridity")
Soil <- c("scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "scalePH")
WaterRetention <- c("scaleSLTPPT", "scaleCLYPPT", "bio10_15_scaled", "scaleAridity", "ScalePET")

groups <- list(
  # ESA = ESA,
  Elevation = Elevation,
  Temperature = Temperature,
  Precip = Precip, 
  Soil = Soil,
  WaterRetention = WaterRetention
)
fgrichness_import_split <- group.importance(fgR_rf, groups) 

######################
# Ordering
#######################
abundance_import_split
abundance_order <- c(1, 3, 4, 6, 2, 5)
richness_import_split
richness_order <- c(4, 5, 6, 1, 2, 3)
biomass_import_split
biomass_order <- c(5,NA, 6,4,3, 2)
fgrichness_import_split
fgrichness_order <- c(NA, 5, 6, 2, 3, 4)

abundance_import_delta <- c(abundance_import_split - max(abundance_import_split, na.rm = TRUE))
richness_import_delta <- c(richness_import_split - max(richness_import_split, na.rm = TRUE))
biomass_import_delta <- c(biomass_import_split - max(biomass_import_split, na.rm = TRUE))
biomass_import_delta <- c(biomass_import_delta[1], NA, biomass_import_delta[2:5])
fgrichness_import_delta <- c(fgrichness_import_split - max(fgrichness_import_split, na.rm = TRUE))
fgrichness_import_delta <- c( NA, fgrichness_import_delta[1:5])

d <- a <- matrix(rep(NA, length = 6*4), nrow = 4, ncol = 6)
colnames(d) <- colnames(a) <- c("ESA","Elevation", "Temperature", "Precipitation", "Soil", "Water Retention")
a[1,] <- richness_order
a[2,] <- abundance_order
a[3,] <- biomass_order
a[4,] <- fgrichness_order
## 4 is most important, 1 is least important
rownames(d) <- rownames(a) <- c("Species Richness", "Abundance", "Biomass", "Functional Richness")


d[1,] <- richness_import_delta
d[2,] <- abundance_import_delta
d[3,] <- biomass_import_delta
d[4,] <- fgrichness_import_delta

d <- round(d, digits = 0)

dat <- melt(a)


dat$X1 <- factor(dat$X1, levels = c( "Functional Richness", "Biomass", "Abundance","Species Richness"))
dat$X2 <- factor(dat$X2, levels = c("ESA","Elevation", "Soil","Precipitation", "Temperature","Water Retention"))

############################

jpeg(file = file.path(figures, "variableImportance_splitGroups.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
p <- VariableImportancePlot(dat, lowColour = "#BCBDDC", highColour = "#25004b", yLab = "Main Models", deltas = d)
p <- p + annotate("text", x = c(1,2,3,4,5, 6), y=1, label = c(d[4,'ESA'],d[4, 'Elevation'],  d[4,'Soil'], d[4,'Precipitation'], d[4,'Temperature'], d[4,'Water Retention'])) +
annotate("text", x = c(1,2,3,4,5,6), y=2, label = c(d[3,'ESA'], d[3, 'Elevation'], d[3,'Soil'], d[3,'Precipitation'], d[3,'Temperature'], d[3,'Water Retention'])) +
annotate("text", x = c(1,2,3,4,5,6), y=3, label = c(d[2,'ESA'], d[2, 'Elevation'], d[2,'Soil'], d[2,'Precipitation'], d[2,'Temperature'], d[2,'Water Retention'])) +
annotate("text", x = c(1,2,3,4,5,6), y=4, label = c(d[1,'ESA'], d[1, 'Elevation'], d[1,'Soil'], d[1,'Precipitation'], d[1,'Temperature'], d[1,'Water Retention']))
p
dev.off()

##### Alternative plot

delta <- melt(d)

alldat <- cbind(dat, delta$value)
names(alldat)[names(alldat) == "delta$value"] <- "delta"

ord <- c(1, 2, 4, 5, 3, 6)

spR <- alldat[alldat$X1 == "Species Richness",]
spR$y <- 4
spR$x <- ord
abund <- alldat[alldat$X1 == "Abundance",]
abund$y <- 3
abund$x <- ord
bmass<- alldat[alldat$X1 == "Biomass",]
bmass$y <- 2
bmass$x <- ord
frichness <- alldat[alldat$X1 == "Functional Richness",]
frichness$y <- 1
frichness$x <- ord



spR <- circleSize(spR)
bmass<- circleSize(bmass)
abund<- circleSize(abund)
frichness<- circleSize(frichness)

all_dat <- rbind(spR, abund, bmass, frichness)

labs <- c("ESA","Elevation","Soil","Precipitation","Temperature","Water\nRetention")

jpeg(file = file.path(figures, "variableImportance_splitGroups_circles.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
par(mar = c(3, 9.5, 1, 5))
plot(-1e+05, -1e+05, ylim = c(0, 5), xlim = c(0.5, 6.5),  
     ylab = "", xlab = "",  xaxt='n', axes = FALSE)
axis(side = 2, cex.axis = 1, labels = levels(all_dat$X1), 
     at = c(1:4), las = 2)
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
fg <- importance(fgR_rf)[,1]

allVars <- unique(c(names(rich), names(abund),names(bio),names(fg)))
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
rich$y <- 4
abund <- getPlottingDf(abund, allVars)
abund$y <- 3
bio <- getPlottingDf(bio, allVars)
bio$y <- 2
fg <- getPlottingDf(fg, allVars)
fg$y <- 1

allDat <- rbind(rich, abund, bio, fg)

labs <- c( "Habitat Cover","Elevation",  "CEC"    ,
          "Clay","Org C","PH"        ,
          "Silt","bio10_7" , "PET"       ,
          "bio10_15", "Snow Months","Aridity Index"   ,
          "bio10_12")

jpeg(file = file.path(figures, "variableImportanceMSE_splitGroups_circles.jpg"), quality = 100, res = 200, width = 2000, height = 1000)

par(mar = c(10, 7, 1, 1))
plot(-1e+05, -1e+05, ylim = c(0, 5), xlim = c(0, length(allVars)),
     ylab = "", xlab = "",  xaxt='n', axes = FALSE)
axis(side = 2, cex.axis = 1, labels = c("Functional Richness", "Biomass", "Abundance", "Species Richness"), 
     at = c(1:4), las = 2)
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
