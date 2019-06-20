########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
library(lme4)
library(car)
library(DHARMa)
library(reshape)
library(Hmisc)
library(ade4)
source(file.path("Functions", "FormatData.R"))
source(file.path("Functions", "Plots.R"))
source(file.path("Functions", "ColourPicker.R"))


#################################################
# 2. Loading in variables
#################################################
models <- "Models"
data_in <-"18_Data"



if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"



files <- list.files(file.path(data_in))
files <- files[grep("Abundance", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

abundance <- read.csv(file.path(data_in, loadin))
abundance <- SiteLevels(abundance) ## 
abundance <- droplevels(abundance)

files <- list.files(file.path(data_in))
files <- files[grep("Biomass", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

biomass <- read.csv(file.path(data_in, loadin))
biomass <- SiteLevels(biomass) ## relevels all land use/habitat variables
biomass <- droplevels(biomass)


files <- list.files(file.path(data_in))
files <- files[grep("Richness", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

richness <- read.csv(file.path(data_in, loadin))
richness <- SiteLevels(richness) ## 
richness <- droplevels(richness)

#############################################
## models
############################################
load(file.path(models, "abundancemodel_functionalgroups.rds"))
load(file.path(models, "biomassmodel_functionalgroups_revised.rds"))
load(file.path(models, "richnessmodel_functionalgroups.rds"))


#############################################
## BIOMASS
############################################
biomassCols <- ColourPicker(biomass$ESA)


newdata <- createNewdata(model = biomass_model, modelFixedEffects = c("ESA","variable", "ScaleElevation",  "scalePH", "scaleCLYPPT", "scaleSLTPPT",  
                                                                "scaleCECSOL", "bio10_4_scaled" , "bio10_12_scaled", "bio10_15_scaled", "SnowMonths_cat", "ScalePETSD"),
                         mainEffect = c("ESA", "variable"), data = biomass)


## Only reference level for snow
newdata <- newdata[newdata$SnowMonths_cat == 0,]

## Predict response and variance
newdata <- predictValues(model = biomass_model, newdata, responseVar = "logValue", 
                         seMultiplier = 1.96, re.form = NA)

#### Get rid of levels not represented by data
fgtokeep <- c("Epi_biomass","Endo_biomass","Ane_biomass")
newdata <- newdata[newdata$variable %in% fgtokeep,]


###############################################
## TRIANGULAR PLOT
##############################################

labelsESA <- as.factor(c("Broadleaf deciduous forest",
              "Needleleaf evergreen forest",
               "Mixed forest","Herbaceous with spare tree/shrub",
               "Shrub","Herbaceous","Production - Herbaceous","Production - Plantation") )

Cols <- biomassCols
# Cols <- ColourPicker(labelsESA)
# This is not in the right order!!

 df <- data.frame(epigeic = newdata$logValue[grep("Epi", newdata$variable)], 
                  endogeics = newdata$logValue[grep("Endo", newdata$variable)],
                  anecics = newdata$logValue[grep("Ane", newdata$variable)])


df$col <- paste0("#", Cols)
row.names(df) <- labelsESA

df[,1:3] <- exp(df[,1:3]) - 1

# Some less than zero
df[which(df[,1] < 0), 1] <- 0
df[which(df[,2] < 0), 2] <- 0
df[which(df[,3] < 0), 3] <- 0


df$total <- rowSums(df[,1:3])



jpeg(file = file.path(figures, "BiomassFGTriangle.jpeg"), quality = 100, res = 200, width = 2000, height = 1000)

t <- triangle.plot(df[,1:3], show.position = FALSE, min3 = c(0, 0, 0), max3 = c(1, 1, 1), cpoint = 0)
points(t, col = df$col, cex = df$total/5 + 0.5, pch = 19)
# Legend for point size
pts <- seq(floor(min(df[,5])), ceiling(max(df[,5])), length.out = 4)
pts <- round(pts)


legend(-1.5,1, legend= pts, pt.cex = pts/5 + 0.5, cex = 0.8, 
       bty = "n", pch = rep(19, times = 4), y.intersp = 2.5, x.intersp =3)
mtext(expression("Total Biomass (grams per" ~ m^{2} ~ ")"), 3, 1.2, adj=0.08)

# legend for colours
legend(0.7,0.8, legend=labelsESA, fill = df$col, cex = 0.8, bty = "n")
dev.off()





##########################################
# Abundance
##########################################

abundanceCols <- ColourPicker(abundance$ESA)


newdata <- createNewdata(model = abundance_model, 
                         modelFixedEffects = c("ESA", "variable", "scalePH", "scaleCLYPPT" ,
                                               "scaleSLTPPT",  "scaleORCDRC", "scaleCECSOL", "bio10_7_scaled", 
                                               "bio10_15_scaled", "SnowMonths_cat","scaleAridity", "ScalePET"),
                         mainEffect = c("ESA", "variable"), data = abundance)


## Only reference level for snow
newdata <- newdata[newdata$SnowMonths_cat == 0,]

## Predict response and variance
newdata <- predictValues(model = abundance_model, newdata, responseVar = "logValue", 
                         seMultiplier = 1.96, re.form = NA)


#### Get rid of levels not represented by data
fgtokeep <- c("Epi_abundance","Endo_abundance","Ane_abundance")
newdata <- newdata[newdata$variable %in% fgtokeep,]

###############################################
## TRIANGULAR PLOT
##############################################

labelsESA <- as.factor(c("Broadleaf deciduous forest", "Broadleaf evergreen forest",
                         "Needleleaf evergreen forest",
                         "Mixed forest","Herbaceous with spare tree/shrub",
                         "Shrub","Herbaceous","Production - Herbaceous","Production - Plantation") )

Cols <- abundanceCols
# Cols <- ColourPicker(labelsESA)
# This is not in the right order!!

df <- data.frame(epigeic = newdata$logValue[grep("Epi", newdata$variable)], 
                 endogeics = newdata$logValue[grep("Endo", newdata$variable)],
                 anecics = newdata$logValue[grep("Ane", newdata$variable)])


df$col <- paste0("#", Cols)
row.names(df) <- labelsESA

df[,1:3] <- exp(df[,1:3]) - 1

# Some less than zero
df[which(df[,1] < 0), 1] <- 0
df[which(df[,2] < 0), 2] <- 0
df[which(df[,3] < 0), 3] <- 0

df$total <- rowSums(df[,1:3])

jpeg(file = file.path(figures, "AbundanceFGTriangle.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
par(mar=c(1,1,1,1))
t <- triangle.plot(df[,1:3], show.position = FALSE, min3 = c(0, 0, 0), max3 = c(1, 1, 1))
points(t, col = df$col, cex = df$total/5, pch = 19)
# legend for colours
legend(0.7,0.8, legend=labelsESA, fill = df$col, cex = 0.8, bty = "n")

# Legend for point size
pts <- seq(floor(min(df[,5])), ceiling(max(df[,5])), length.out = 4)
pts <- round(pts)
legend(-1.5,0.8, legend= pts, fill = "black",pt.cex = pts/5, cex = 0.8, 
       bty = "n", pch = rep(19, times = 4), y.intersp = 2.7, x.intersp =3)
mtext("Total Abundance (ind./m2)", 3, -3, adj=0.1)

dev.off()


##########################################
# Species Richness
##########################################

richnessCols <- ColourPicker(richness$ESA)


newdata <- createNewdata(model = richness_model, 
                         modelFixedEffects = c("ESA", "variable", "scalePH", 
                                               "scaleSLTPPT", 
                                               "bio10_15_scaled", "SnowMonths_cat", "ScalePET"),
                         mainEffect = c("ESA", "variable"), data = abundance)


## Only reference level for snow
newdata <- newdata[newdata$SnowMonths_cat == 0,]

## Predict response and variance
newdata <- predictValues(model = richness_model, newdata, responseVar = "value", 
                         seMultiplier = 1.96, re.form = NA)


fgtokeep <- c("Epi_richness","Endo_richness","Ane_richness")
newdata <- newdata[newdata$variable %in% fgtokeep,]

###############################################
## TRIANGULAR PLOT
##############################################

labelsESA <- as.factor(c("Broadleaf deciduous forest", "Broadleaf evergreen forest",
                         "Needleleaf evergreen forest",
                         "Mixed forest","Herbaceous with spare tree/shrub",
                         "Shrub","Herbaceous","Production - Herbaceous","Production - Plantation") )

Cols <- richnessCols
# Cols <- ColourPicker(labelsESA)
# This is not in the right order!!

df <- data.frame(epigeic = newdata$value[grep("Epi", newdata$variable)], 
                 endogeics = newdata$value[grep("Endo", newdata$variable)],
                 anecics = newdata$value[grep("Ane", newdata$variable)])


df$col <- paste0("#", Cols)
row.names(df) <- labelsESA

df[,1:3] <- exp(df[,1:3])

df$total <- rowSums(df[,1:3])

jpeg(file = file.path(figures, "RichnessFGTriangle.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
par(mar=c(1,1,1,1))
t <- triangle.plot(df[,1:3], show.position = FALSE, min3 = c(0, 0, 0), max3 = c(1, 1, 1))
points(t, col = df$col, cex = df$total/5, pch = 19)
# legend for colours
legend(0.7,0.8, legend=labelsESA, fill = df$col, cex = 0.8, bty = "n")

# Legend for point size
pts <- seq(floor(min(df[,5])), ceiling(max(df[,5])), length.out = 4)
pts <- round(pts)
legend(-1.5,0.8, legend= pts, fill = "black",pt.cex = pts/5, cex = 0.8, 
       bty = "n", pch = rep(19, times = 4), y.intersp = 2.7, x.intersp =3)
mtext("Total Abundance (ind./m2)", 3, -3, adj=0.1)

dev.off()

