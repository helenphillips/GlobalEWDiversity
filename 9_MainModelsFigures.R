########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
source(file.path("Functions", "Plots.R"))
source(file.path("Functions", "ColourPicker.R"))

library(lme4)
library(Hmisc)
library(maps)
library(maptools)

source("Functions/FormatData.R")

#################################################
# 2. Loading in variables
#################################################

data_in <-"8_Data"

files <- list.files(file.path(data_in))
files <- files[grep("sitesRichness", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
# loadin <- files[grep(date, files)]
richness <- read.csv(file.path(data_in, paste("sitesRichness_",date,".csv", sep = "")))


files <- list.files(file.path(data_in))
files <- files[grep("Biomass", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
# loadin <- files[grep(date, files)]
biomass <- read.csv(file.path(data_in, paste("sitesBiomass_",date,".csv", sep = "")))

files <- list.files(file.path(data_in))
files <- files[grep("Abundance", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
# loadin <- files[grep(date, files)]
abundance <- read.csv(file.path(data_in, paste("sitesAbundance_",date,".csv", sep = "")))


if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

#################################################
# 3. Load in data
#################################################


richness <- droplevels(SiteLevels(richness))
biomass <- droplevels(SiteLevels(biomass))
abundance <- droplevels(SiteLevels(abundance))

#################################################
# 4. Load in models
#################################################

models <- "Models"


load(file.path(models, "richnessmodel.rds"))
load(file.path(models, "biomassmodel_full.rds"))
load(file.path(models, "abundancemodel_full.rds"))

#################################################
# 5. Pick colors for figures
#################################################

biomassCols <- ColourPicker(biomass$ESA)
richnessCols <- ColourPicker(richness$ESA)
abundanceCols <- ColourPicker(abundance$ESA)


############################################
# SPECIES RICHNES
#############################################


jpeg(file = file.path(figures, "richness_ESA.jpg"), quality = 100, res = 200, width = 2000, height = 1300)
plotSingle(model= richness_model, backTransform = TRUE, family = "poisson",
           modelFixedEffs = c("scalePH","scaleCLYPPT", "scaleSLTPPT" , "scaleCECSOL" ,
                                "scaleORCDRC", "bio10_7_scaled" , "bio10_15_scaled" , "SnowMonths_cat",  
                                "scaleAridity", "ScalePET", "ESA", "scaleElevation"),
           Effect1 = "ESA", 
           responseVar = "SpeciesRichness", seMultiplier = 1, data = richness, cols = richnessCols, 
           legend.position, ylabel = "Species Richness", xlabel = "", otherContEffectsFun = "median")
dev.off()

#############################################
# ABUNDANCE
###############################################

jpeg(file = file.path(figures, "abundance_ESA.jpg"), quality = 100, res = 200, height = 1000, width = 1500, pointsize = 12)
plotSingle(model= abundance_model, Effect1 = "ESA",
           modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleSLTPPT" , "scaleCECSOL" ,
                                "scaleORCDRC", "bio10_7_scaled" , "bio10_15_scaled" ,"SnowMonths_cat",  
                                "scaleAridity" , "ScalePET", "ESA" , "ScaleElevation"),
           responseVar = "logAbundance", seMultiplier = 1, data = abundance, cols = abundanceCols, 
           legend.position = NA, ylabel = "log(Abundance)", xlabel = "", otherContEffectsFun = "median")
dev.off()

#################################################

jpeg(file = file.path(figures, "biomass_ESA.jpg"), quality = 100, res = 200, width = 2000, height = 1300)
plotSingle(model= biomass_model, Effect1 = "ESA", 
           modelFixedEffs = c("scalePH" ,"scaleCLYPPT" , "scaleSLTPPT" , "scaleORCDRC",
                                "scaleCECSOL" , "bio10_7_scaled" ,"bio10_12_scaled", "bio10_15_scaled",  
                                "ScalePET" , "SnowMonths_cat", "ESA"),
           responseVar = "logBiomass", seMultiplier = 1, data = biomass, cols = biomassCols, 
           legend.position = NA, ylabel = "log(Biomass)", xlabel = "", otherContEffectsFun = "median")
dev.off()



##############################################
## FOr the manuscript
##############################################
jpeg(file = file.path(figures, "HabitatCover_Richness+Abundance.jpg"), quality = 100, res = 200, width = 1500, height = 2000)

par(mfrow = c(2, 1))
par(mar = c(15, 3, 1, 1))
plotSingle(model= richness_model, backTransform = TRUE, family = "poisson",
           modelFixedEffs = c("scalePH","scaleCLYPPT", "scaleSLTPPT" , "scaleCECSOL" ,
                              "scaleORCDRC", "bio10_7_scaled" , "bio10_15_scaled" , "SnowMonths_cat",  
                              "scaleAridity", "ScalePET", "ESA", "scaleElevation"),
           Effect1 = "ESA", 
           responseVar = "SpeciesRichness", seMultiplier = 1, data = richness, cols = richnessCols, 
           legend.position, ylabel = "Species Richness", xlabel = "", otherContEffectsFun = "median")
mtext("(a)", side = 3, line = 0, at = 0, adj = 0.1)

plotSingle(model= abundance_model, Effect1 = "ESA",
           modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleSLTPPT" , "scaleCECSOL" ,
                              "scaleORCDRC", "bio10_7_scaled" , "bio10_15_scaled" ,"SnowMonths_cat",  
                              "scaleAridity" , "ScalePET", "ESA" , "ScaleElevation"),
           responseVar = "logAbundance", seMultiplier = 1, data = abundance, cols = abundanceCols, 
           legend.position = NA, ylabel = "log(Abundance)", xlabel = "", otherContEffectsFun = "median")
mtext("(b)", side = 3, line = 0, at = 0, adj = 0.1)



dev.off()

##################################################
## MAP
##################################################
all_data <-"7_Data"

files <- list.files(file.path(all_data))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
all_data <- read.csv(file.path(all_data, paste("sites_",date,".csv", sep = "")))


studies1 <- as.vector(unique(richness_model@frame$Study_Name))
studies2 <- as.vector(unique(abundance_model@frame$Study_Name))
studies3 <- as.vector(unique(biomass_model@frame$Study_Name))

all_studies <- c(studies1, studies2, studies3)
all_studies <- unique(all_studies)



all_studies <- all_data[all_data$Study_Name %in% all_studies,]




coord<-aggregate(cbind(all_studies$Longitude__Decimal_Degrees, all_studies$Latitude__decimal_degrees), list(all_studies$Study_Name), mean)
coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")




#pdf(file = file.path(figures, "Map_alldata.pdf"), height = 4)
jpeg(filename = file.path(figures, "Map_modelledData.jpg"), quality = 100, res = 300, width = 2000, height = 2000)
mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)
dev.off()


nrow(all_studies)
length(unique(all_studies$file))
length(unique(all_studies$Study_Name))
unique(all_studies$Country)
