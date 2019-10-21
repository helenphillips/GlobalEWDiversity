########################################################
# 0. Set Working Directory
########################################################
if(Sys.info()["nodename"] == "IDIVNB179"){
   setwd("C:\\USers\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "TSGIS02"){
   setwd("C:/sWorm/EarthwormAnalysis")
}
#################################################
# 1. Loading libraries
#################################################
source(file.path("Functions", "Plots.R"))
source(file.path("Functions", "ColourPicker.R"))
source(file.path("Functions", "FormatData.R"))
source(file.path("Functions", "cornerlabel2.R"))
 
 
library(lme4)
library(Hmisc)
library(maps)
library(maptools)

#######################################
# variables
######################################
 
wide_cm <- 12
wide_inch <- 4.75
point_size <- 7
plotlabcex <- 0.5
resdpi <- 300

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

if(!dir.exists("9_Data")){
  dir.create("9_Data")
}
data_out <- "9_Data"

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


load(file.path(models, "richnessmodel_revised.rds"))
load(file.path(models, "biomassmodel_full_revised.rds"))
load(file.path(models, "abundancemodel_full_revised.rds"))

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
mtext("A", side = 3, line = 0, at = 0, adj = 0.1, font = 2)

plotSingle(model= abundance_model, Effect1 = "ESA",
           modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleSLTPPT" , "scaleCECSOL" ,
                              "scaleORCDRC", "bio10_7_scaled" , "bio10_15_scaled" ,"SnowMonths_cat",  
                              "scaleAridity" , "ScalePET", "ESA" , "ScaleElevation"),
           responseVar = "logAbundance", seMultiplier = 1, data = abundance, cols = abundanceCols, 
           legend.position = NA, ylabel = "(ln-) Abundance", xlabel = "", otherContEffectsFun = "median")
mtext("B", side = 3, line = 0, at = 0, adj = 0.1, font = 2)



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



## Save for future use
write.csv(all_studies, file = file.path(data_out, "sWorm_CompleteDataSet.csv"), row.names = FALSE)
## 

coord<-aggregate(cbind(all_studies$Longitude__Decimal_Degrees, all_studies$Latitude__decimal_degrees), list(all_studies$Study_Name), mean)
coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
#proj4string(dsSPDF)<-CRS("+proj=longlat")
proj4string(dsSPDF)<-CRS("+init=ESRI:54030")


#pdf(file = file.path(figures, "Map_alldata.pdf"), height = 4)
jpeg(filename = file.path(figures, "Map_modelledData.jpg"), quality = 100, res = 300, width = 2000, height = 2000)
mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)
plotrix::corner.label(label = "(a)", x = -1, y = 1, cex = plotlabcex)

dev.off()

#############################################
## Alternate map
##############################################


all(coord$X %in% names(table(all_studies$Study_Name)))

studyN <- data.frame(table(all_studies$Study_Name))

coord <- merge(coord, studyN, all.x = TRUE, by.x = "X", by.y = "Var1")


# coord$X<-coord$Group.1
# coord<-coord[2:4]
# coord$Group.1 <- NULL
names(coord)<-c("X", "Long", "Lat", "nSites")


coord$size <- ((coord$nSites-min(coord$nSites))/(max(coord$nSites)-min(coord$nSites)) * 2) + 0.5


dsSPDF<-SpatialPointsDataFrame(coord[,2:3], data.frame(coord[,1:5]))
proj4string(dsSPDF)<-CRS("+proj=longlat")

transpBlack <- rgb(0, 0, 0, alpha = 0.4, names = NULL, maxColorValue = 1)

# jpeg(filename = file.path(figures, "Map_modelledData_nsites.jpg"), quality = 100, res = 300, width = 2000, height = 2000)
pdf(file.path(figures, "Map_modelledData_nsites.pdf"),width= wide_inch, height= wide_inch/2, pointsize = point_size)

mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col=transpBlack, bg = transpBlack, cex= coord$size, pch=19)
corner.label2(label = "A", x = -1, y = 1, cex = plotlabcex, font = 2)


sizes <- c(1, 50, 100, 150, 200, 250)
cexsizes <- ((sizes-min(coord$nSites))/(max(coord$nSites)-min(coord$nSites)) * 2) + 0.5

legend(x = -170, y = 2, legend = sizes, pch = 19, pt.cex =cexsizes, bty="n", cex = 0.7, 
       y.intersp = c(1, 1, 1, 1.05, 1.1, 1.18),
       x.intersp = c(1.19),
       title = "Number Of Sites")
dev.off()



########################################
nrow(all_studies)
length(unique(all_studies$file))
length(unique(all_studies$Study_Name))
unique(all_studies$Country)
