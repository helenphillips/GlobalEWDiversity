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
#############################################




plotInteraction(model = richness_model, Effect1 = "scalePH", Effect2 = "scaleORCDRC",
                modelFixedEffs = c("scalePH", "scaleORCDRC", "bio10_1_scaled", "bio10_4_scaled", 
                                   "bio10_12_scaled", "bio10_15_scaled", "ESA", "scaleCECSOL"),
                responseVar = "SpeciesRichness", seMultiplier = 1.96,
                data = richness, cols = "black", legend.position = "topleft",
                ylabel = "", xlabel = "")



#############################################
# ABUNDANCE
###############################################

# pdf(file = file.path(figures, "Abundance_ph.pdf"), height = 4)

jpeg(file = file.path(figures, "Abundance_ESA.jpg"), quality = 100, res = 200, height = 1000, width = 1500, pointsize = 12)
plotSingle(model= abundance_model, modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
            "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
           Effect1 = "ESA", 
           responseVar = "logAbundance", seMultiplier = 1, data = abundance, cols = abundanceCols, 
           legend.position, ylabel = "log(Abundance)", xlabel = "", otherContEffectsFun = "median")
dev.off()

plotInteraction(model = abundance_model, Effect1 = "scalePH", Effect2 = "scaleCECSOL",
                modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleCECSOL", "scaleORCDRC", 
                                   "ESA", "bio10_1_scaled", "bio10_4_scaled", "bio10_12_scaled", "bio10_15_scaled"),
                responseVar = "logAbundance", seMultiplier = 1.96,
                data = abundance, cols = "black", legend.position = "topleft",
                ylabel = "", xlabel = "")


#################################################
# 6. Figures
#################################################
### Interactions
#pdf(file = file.path(figures, "Biomass_ph.pdf"), height = 4)
jpeg(file = file.path(figures, "Biomass_ph.jpg"), quality = 100, res = 200, height = 1000, width = 2000)
plotSingle(model= biomass_model, modelFixedEffs = c("scalePH", "ESA", "bio10_14_scaled", "bio10_13_scaled", "bio10_5_scaled"), 
           Effect1 = "scalePH", 
                responseVar = "logBiomass", seMultiplier = 1, data = biomass, cols = "000000", 
                legend.position, ylabel = "log(Wet Biomass)", xlabel = "pH", otherContEffectsFun = "median")
dev.off()  


# plotInteraction(model = b2a, Effect1 = "scalePH", Effect2 = "LU_Mgmt", 
#                 modelFixedEffs = c("intensity", "LU_Mgmt", "scalePH"),
#                 responseVar = "logBiomass", seMultiplier = 1.96, 
#                 data = biomass, cols = luMgmtCols, legend.position = "topleft", 
#                 ylabel = "", xlabel = "")
  



# pdf(file = file.path(figures, ""), height = 4)
jpeg(file = file.path(figures, "Biomass_ESA.jpg"), quality = 100, res = 200, height = 1000, width = 1500, pointsize = 12)
plotSingle(model= b4a, modelFixedEffs = c("scalePH", "ESA", "bio10_14_scaled", "bio10_13_scaled", "bio10_5_scaled"),
           Effect1 = "ESA", 
           responseVar = "logBiomass", seMultiplier = 1, data = biomass, cols = habitCols, 
           legend.position, ylabel = "log(Wet Biomass)", xlabel = "", otherContEffectsFun = "median")
dev.off()


## Climate
jpeg(file = file.path(figures, "Biomass_bio10_14.jpg"), quality = 100, res = 200, height = 1000, width = 1500, pointsize = 12)

plotSingle(model= b4a, modelFixedEffs = c("scalePH", "ESA", "bio10_14_scaled", "bio10_13_scaled", "bio10_5_scaled"),
           Effect1 = "bio10_14_scaled", 
           responseVar = "logBiomass", seMultiplier = 1, data = biomass, cols = "000000", 
           legend.position, ylabel = "log(Wet Biomass)", xlabel = "Precip. of Driest Month", otherContEffectsFun = "median")
dev.off()

jpeg(file = file.path(figures, "Biomass_bio10_13.jpg"), quality = 100, res = 200, height = 1000, width = 1500, pointsize = 12)
plotSingle(model= b4a, modelFixedEffs = c("scalePH", "ESA", "bio10_14_scaled", "bio10_13_scaled", "bio10_5_scaled"),
           Effect1 = "bio10_13_scaled", 
           responseVar = "logBiomass", seMultiplier = 1, data = biomass, cols = "000000", 
           legend.position, ylabel = "log(Wet Biomass)", xlabel = "Precip. of Wettest Month", otherContEffectsFun = "median")
dev.off()


jpeg(file = file.path(figures, "Biomass_bio10_5.jpg"), quality = 100, res = 200, height = 1000, width = 1500, pointsize = 12)
plotSingle(model= b4a, modelFixedEffs = c("scalePH", "ESA", "bio10_14_scaled", "bio10_13_scaled", "bio10_5_scaled"),
           Effect1 = "bio10_5_scaled", 
           responseVar = "logBiomass", seMultiplier = 1, data = biomass, cols = "000000", 
           legend.position, ylabel = "log(Wet Biomass)", xlabel = "Max temp. of warmest month", otherContEffectsFun = "median")
dev.off()


#################################################3
## MAP
#################################################3

studies <- unique(a2a@frame$Study_Name)
abundance2 <- abundance[abundance$Study_Name %in% studies,]


coord<-aggregate(cbind(abundance2$Longitude__Decimal_Degrees, abundance2$Latitude__decimal_degrees), list(abundance2$Study_Name), mean)
## Five don't have coordinates yet
coord <- coord[complete.cases(coord),]


coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")



studiesb <- unique(b2a@frame$Study_Name)



#pdf(file = file.path(figures, "Map_alldata.pdf"), height = 4)
jpeg(filename = file.path(figures, "Map_abundance.jpg"), quality = 100, res = 300, width = 2000, height = 2000)
mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)
dev.off()
