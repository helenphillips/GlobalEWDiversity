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

#################################################
# 2. Variables
#################################################


data_in <- "4_Data"
models <- "Models"
figures <- "Figures"
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

load(file.path(models, "richnessmodel_full.rds"))

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
files <- files[grep("sites+FGRichness_", files)]
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
                      "scaleORCDRC" , "bio10_1_scaled" , "bio10_15_scaled" , "SnowMonths_cat" ,
                      "scaleAridity" , "ScalePETSD", "ESA")

abundance_rf <- randomForest(y = abundance$logAbundance, x = abundance[,names(abundance) %in% AbundmainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)

ESA <- "ESA"
Temperature <- c("bio10_1_scaled", "ScalePETSD")
Precip <- c("bio10_15_scaled" , "SnowMonths_cat", "scaleAridity")
Soil <- c("scalePH", "scaleCLYPPT", "scaleCECSOL", 
          "scaleSLTPPT", "scaleORCDRC")
WaterRetention <- c("scaleCLYPPT", "bio10_15_scaled", "ScalePETSD", "scaleAridity")

groups <- list(
  ESA = ESA,
  Temperature = Temperature,
  Precip = Precip,
  Soil = Soil,
  WaterRetention = WaterRetention
)
abundance_import_split <- group.importance(abundance_rf, groups) 

# RICHNESS
richness_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL",  
                          "scaleORCDRC", "bio10_4_scaled", "bio10_15_scaled", "SnowMonths_cat",  
                          "scaleAridity", "ScalePET", "ESA")

spR_rf <- randomForest(y = richness$SpeciesRichness, x = richness[,names(richness) %in% richness_mainEffects], 
                       ntree=501, importance=TRUE, proximity = TRUE)



ESA <- "ESA"
Temperature <- c("bio10_4_scaled", "ScalePET")
Precip <- c("bio10_15_scaled", "SnowMonths_cat", "scaleAridity")
Soil <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleCECSOL", "scaleORCDRC")
WaterRetention <- c("scaleCLYPPT","scaleSLTPPT", "bio10_15_scaled", "scaleAridity", "ScalePET")

groups <- list(
  ESA = ESA,
  Temperature = Temperature,
  Precip = Precip, 
  Soil = Soil,
  WaterRetention = WaterRetention
)
richness_import_split <- group.importance(spR_rf, groups) 

# BIOMASS
biomass_mainEffects <- c("scalePH", "scaleCLYPPT", "scaleSLTPPT", "scaleORCDRC", "scaleCECSOL", "bio10_12_scaled",  
                         "bio10_15_scaled", "SnowMonths_cat", "ESA")

bioM_rf <- randomForest(y = biomass$logBiomass, x = biomass[,names(biomass) %in% biomass_mainEffects], 
                        ntree=501, importance=TRUE, proximity = TRUE)

ESA <- "ESA"
Precip <- c("bio10_12_scaled", "bio10_15_scaled", "SnowMonths_cat")
Soil <- c("scalePH","scaleORCDRC", "scaleCECSOL","scaleCLYPPT", "scaleSLTPPT")
WaterRetention <- c("scaleCLYPPT", "scaleSLTPPT", "bio10_12_scaled", "bio10_15_scaled")

groups <- list(
  ESA = ESA,
  Precip = Precip,
  Soil = Soil,
  WaterRetention = WaterRetention
)
biomass_import_split <- group.importance(bioM_rf, groups) 


## Functional Richness

# RICHNESS
fgrichness_mainEffects <- c("scaleSLTPPT", "scaleORCDRC", "bio10_4_scaled", "bio10_15_scaled",  
                              "SnowMonths_cat", "scaleAridity", "ScalePET")

fgR_rf <- randomForest(y = fgrichness$FGRichness, x = fgrichness[,names(fgrichness) %in% fgrichness_mainEffects], 
                       ntree=501, importance=TRUE, proximity = TRUE)



# ESA <- "ESA"
Temperature <- c("bio10_4_scaled", "ScalePET")
Precip <- c("bio10_15_scaled", "SnowMonths_cat", "scaleAridity")
Soil <- c("scaleSLTPPT", "scaleORCDRC")
WaterRetention <- c("scaleSLTPPT", "bio10_15_scaled", "scaleAridity", "ScalePET")

groups <- list(
  # ESA = ESA,
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
abundance_order <- c(1, 3, 5, 2, 4)
richness_import_split
richness_order <- c(4, 5, 2, 1, 3)
biomass_import_split
biomass_order <- c(4,NA, 5,2,3)
fgrichness_import_split
fgrichness_order <- c(NA, 5, 3, 2, 4)

abundance_import_delta <- c(abundance_import_split - max(abundance_import_split, na.rm = TRUE))
richness_import_delta <- c(richness_import_split - max(richness_import_split, na.rm = TRUE))
biomass_import_delta <- c(biomass_import_split - max(biomass_import_split, na.rm = TRUE))
biomass_import_delta <- c(biomass_import_delta[1], NA, biomass_import_delta[2:4])
fgrichness_import_delta <- c(fgrichness_import_split - max(fgrichness_import_split, na.rm = TRUE))
fgrichness_import_delta <- c( NA, fgrichness_import_delta[1:4])

d <- a <- matrix(rep(NA, length = 5*4), nrow = 4, ncol = 5)
colnames(d) <- colnames(a) <- c("ESA", "Temperature", "Precipitation", "Soil", "Water Retention")
a[1,] <- richness_order
a[2,] <- abundance_order
a[3,] <- biomass_order
a[4,] <- fgrichness_order
## 4 is most important, 1 is least important
rownames(d) <- rownames(a) <- c("SpeciesRichness", "Abundance", "Biomass", "Functional Richness")


d[1,] <- richness_import_delta
d[2,] <- abundance_import_delta
d[3,] <- biomass_import_delta
d[4,] <- fgrichness_import_delta

dat <- melt(a)

dat$X1 <- factor(dat$X1, levels = c( "Functional Richness", "Biomass", "Abundance","SpeciesRichness"))
dat$X2 <- factor(dat$X2, levels = c("ESA","Soil","Precipitation", "Temperature","Water Retention"))

############################

jpeg(file = file.path(figures, "variableImportance_splitGroups.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
p <- VariableImportancePlot(dat, lowColour = "#BCBDDC", highColour = "#25004b", yLab = "Main Models", deltas = d)
p <- p + annotate("text", x = c(1,2,3,4,5), y=1, label = c(d[4,'ESA'], d[4,'Soil'], d[4,'Precipitation'], d[4,'Temperature'], d[4,'Water Retention'])) +
annotate("text", x = c(1,2,3,4,5), y=2, label = c(d[3,'ESA'], d[3,'Soil'], d[3,'Precipitation'], d[3,'Temperature'], d[3,'Water Retention'])) +
annotate("text", x = c(1,2,3,4,5), y=3, label = c(d[2,'ESA'], d[2,'Soil'], d[2,'Precipitation'], d[2,'Temperature'], d[2,'Water Retention'])) +
annotate("text", x = c(1,2,3,4,5), y=4, label = c(d[1,'ESA'], d[1,'Soil'], d[1,'Precipitation'], d[1,'Temperature'], d[1,'Water Retention']))
p
dev.off()

##### Alternative plot

delta <- melt(d)

alldat <- cbind(dat, delta$value)
names(alldat)[names(alldat) == "delta$value"] <- "delta"

spR <- alldat[alldat$X1 == "SpeciesRichness",]
spR$y <- 4
spR$x <- c(1, 3, 4, 2, 5)
abund <- alldat[alldat$X1 == "Abundance",]
abund$y <- 3
abund$x <- c(1, 3, 4, 2, 5)
bmass<- alldat[alldat$X1 == "Biomass",]
bmass$y <- 2
bmass$x <- c(1, 3, 4, 2, 5)
frichness <- alldat[alldat$X1 == "Functional Richness",]
frichness$y <- 1
frichness$x <- c(1, 3, 4, 2, 5)

circleSize <- function(dat){
  dat$size <- 9
  onepercent <- min(dat$delta, na.rm = TRUE) / 100
  dat$percentChange <- dat$delta / onepercent
  dat$size <- dat$size * ((100 - dat$percentChange) / 100)
  dat$size <- dat$size + 1
  return(dat)
}

spR <- circleSize(spR)
bmass<- circleSize(bmass)
abund<- circleSize(abund)
frichness<- circleSize(frichness)

all_dat <- rbind(spR, abund, bmass, frichness)

jpeg(file = file.path(figures, "variableImportance_splitGroups_circles.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
par(mar = c(3, 9.5, 1, 3))
plot(-1e+05, -1e+05, ylim = c(0, 5), xlim = c(0.5, 5.5),  
     ylab = "", xlab = "",  xaxt='n', axes = FALSE)
Axis(side = 2, cex.axis = 1, labels = levels(all_dat$X1), 
     at = c(1:4), las = 2)
points(all_dat$x, all_dat$y, pch = 19, cex = all_dat$size, ylim = c(0, 5))
axis(side=1, at = 1:5, labels = levels(all_dat$X2), las=1, cex.axis = 1)
mtext("Model", side = 2, line = 6, cex = 2)
dev.off()
