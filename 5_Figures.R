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

data_in <-"4_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
# loadin <- files[grep(date, files)]

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

#################################################
# 3. Load in data
#################################################

richness <- read.csv(file.path(data_in, paste("sitesRichness_",date,".csv", sep = "")))
biomass <- read.csv(file.path(data_in, paste("sitesBiomass_",date,".csv", sep = "")))
abundance <- read.csv(file.path(data_in, paste("sitesAbundance_",date,".csv", sep = "")))

richness <- droplevels(SiteLevels(richness))
biomass <- droplevels(SiteLevels(biomass))
abundance <- droplevels(SiteLevels(abundance))



#################################################
# 4. Load in models
#################################################

models <- "Models"


load(file.path(models, "richnessmodel_full.rds"))
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
####(#########################################
jpeg(file = file.path(figures, "richness_pHOrgc.jpg"))
plotInteraction(model = richness_model, Effect1 = "scalePH", Effect2 = "scaleORCDRC",
                modelFixedEffs = c("scalePH", "scaleORCDRC", "bio10_1_scaled", "bio10_4_scaled", 
                                   "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
                responseVar = "SpeciesRichness", seMultiplier = 1.96,
                data = richness, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "richness_Bio12Bio15.jpg"))
plotInteraction(model = richness_model, Effect1 = "bio10_12_scaled", Effect2 = "bio10_15_scaled",
                modelFixedEffs = c("scalePH", "scaleORCDRC", "bio10_1_scaled", "bio10_4_scaled", 
                                   "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
                responseVar = "SpeciesRichness", seMultiplier = 1.96,
                data = richness, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "richness_Bio1Bio4.jpg"))
plotInteraction(model = richness_model, Effect1 = "bio10_1_scaled", Effect2 = "bio10_4_scaled",
                modelFixedEffs = c("scalePH", "scaleORCDRC", "bio10_1_scaled", "bio10_4_scaled", 
                                   "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
                responseVar = "SpeciesRichness", seMultiplier = 1.96,
                data = richness, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "richness_ESA.jpg"))
plotSingle(model= richness_model, 
           modelFixedEffs = c("scalePH", "scaleORCDRC", "bio10_1_scaled", "bio10_4_scaled", 
                              "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
           Effect1 = "ESA", 
           responseVar = "SpeciesRichness", seMultiplier = 1, data = abundance, cols = abundanceCols, 
           legend.position, ylabel = "Species Richness", xlabel = "", otherContEffectsFun = "median")
dev.off()

jpeg(file = file.path(figures, "richness_CEC.jpg"))
plotSingle(model= richness_model, 
           modelFixedEffs = c("scalePH", "scaleORCDRC", "bio10_1_scaled", "bio10_4_scaled", 
                              "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
           Effect1 = "scaleCECSOL", 
           responseVar = "SpeciesRichness", seMultiplier = 1, data = abundance, cols = abundanceCols, 
           legend.position, ylabel = "Species Richness", xlabel = "CEC", otherContEffectsFun = "median")
dev.off()
#############################################
# ABUNDANCE
###############################################

# pdf(file = file.path(figures, "Abundance_ph.pdf"), height = 4)
jpeg(file = file.path(figures, "abundance_Bio12Bio15.jpg"))
plotInteraction(model = abundance_model, Effect1 = "bio10_12_scaled", Effect2 = "bio10_15_scaled",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "abundance_Bio4Bio12.jpg"))
plotInteraction(model = abundance_model, Effect1 = "bio10_4_scaled", Effect2 = "bio10_12_scaled",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "abundance_Bio1Bio15.jpg"))
plotInteraction(model = abundance_model, Effect1 = "bio10_1_scaled", Effect2 = "bio10_15_scaled",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "abundance_Bio1Bio12.jpg"))
plotInteraction(model = abundance_model, Effect1 = "bio10_1_scaled", Effect2 = "bio10_12_scaled",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "abundance_Bio1Bio4.jpg"))
plotInteraction(model = abundance_model, Effect1 = "bio10_1_scaled", Effect2 = "bio10_4_scaled",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "abundance_CECOrgC.jpg"))
plotInteraction(model = abundance_model, Effect1 = "scaleCECSOL", Effect2 = "scaleORCDRC",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "abundance_ClayCEC.jpg"))
plotInteraction(model = abundance_model, Effect1 = "scaleCLYPPT", Effect2 = "scaleCECSOL",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "abundance_pHOrgC.jpg"))
plotInteraction(model = abundance_model, Effect1 = "scalePH", Effect2 = "scaleORCDRC",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "abundance_pHCEC.jpg"))
plotInteraction(model = abundance_model, Effect1 = "scalePH", Effect2 = "scaleCECSOL",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "abundance_ESA.jpg"), quality = 100, res = 200, height = 1000, width = 1500, pointsize = 12)
plotSingle(model= abundance_model, modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
            "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
           Effect1 = "ESA", 
           responseVar = "logAbundance", seMultiplier = 1, data = abundance, cols = abundanceCols, 
           legend.position, ylabel = "log(Abundance)", xlabel = "", otherContEffectsFun = "median")
dev.off()




#################################################
# BIOMASS
#################################################
jpeg(file = file.path(figures, "biomass_Bio4Bio12.jpg"))
plotInteraction(model = biomass_model, Effect1 = "bio10_4_scaled", Effect2 = "bio10_12_scaled",
                modelFixedEffs = c("scalePH", "scaleORCDRC", "scaleCLYPPT", "bio10_1_scaled", "bio10_4_scaled", 
                                   "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
                responseVar = "logBiomass", seMultiplier = 1.96,
                data = biomass, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "biomass_Bio1Bio12.jpg"))
plotInteraction(model = biomass_model, Effect1 = "bio10_1_scaled", Effect2 = "bio10_12_scaled",
                modelFixedEffs = c("scalePH", "scaleORCDRC", "scaleCLYPPT", "bio10_1_scaled", "bio10_4_scaled", 
                                   "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
                responseVar = "logBiomass", seMultiplier = 1.96,
                data = biomass, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "biomass_ClayOrgC.jpg"))
plotInteraction(model = biomass_model, Effect1 = "scaleCLYPPT", Effect2 = "scaleORCDRC",
                modelFixedEffs = c("scalePH", "scaleORCDRC", "scaleCLYPPT", "bio10_1_scaled", "bio10_4_scaled", 
                                   "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
                responseVar = "logBiomass", seMultiplier = 1.96,
                data = biomass, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "biomass_pHCEC.jpg"))
plotInteraction(model = biomass_model, Effect1 = "scalePH", Effect2 = "scaleCECSOL",
                modelFixedEffs = c("scalePH", "scaleORCDRC", "scaleCLYPPT", "bio10_1_scaled", "bio10_4_scaled", 
                                   "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
                responseVar = "logBiomass", seMultiplier = 1.96,
                data = biomass, cols = c("white", "red"), legend.position = "topleft",
                ylabel = "", xlabel = "")
dev.off()

jpeg(file = file.path(figures, "biomass_ESA.jpg"))
plotSingle(model= biomass_model, 
           modelFixedEffs = c("scalePH", "scaleORCDRC", "scaleCLYPPT", "bio10_1_scaled", "bio10_4_scaled", 
                              "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
           Effect1 = "ESA", 
           responseVar = "logBiomass", seMultiplier = 1, data = biomass, cols = biomassCols, 
           legend.position, ylabel = "log(Biomass)", xlabel = "", otherContEffectsFun = "median")
dev.off()

##################################################
## MAP
##################################################
all_data <-"3_Data"

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
