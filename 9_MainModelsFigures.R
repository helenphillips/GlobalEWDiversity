########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
source(file.path("Functions", "Plots.R"))
source(file.path("Functions", "ColourPicker.R"))
source(file.path("Functions", "FormatData.R"))

library(lme4)
library(Hmisc)
library(maps)
library(maptools)


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
mtext("(a)", side = 3, line = 0, at = 0, adj = 0.1)

plotSingle(model= abundance_model, Effect1 = "ESA",
           modelFixedEffs = c("scalePH", "scaleCLYPPT", "scaleSLTPPT" , "scaleCECSOL" ,
                              "scaleORCDRC", "bio10_7_scaled" , "bio10_15_scaled" ,"SnowMonths_cat",  
                              "scaleAridity" , "ScalePET", "ESA" , "ScaleElevation"),
           responseVar = "logAbundance", seMultiplier = 1, data = abundance, cols = abundanceCols, 
           legend.position = NA, ylabel = "(ln-) Abundance", xlabel = "", otherContEffectsFun = "median")
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

jpeg(filename = file.path(figures, "Map_modelledData_nsites.jpg"), quality = 100, res = 300, width = 2000, height = 2000)
mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col=transpBlack, bg = transpBlack, cex= coord$size, pch=19)
plotrix::corner.label(label = "(a)", x = -1, y = 1, cex = plotlabcex)


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


####################
# FROM JO
####################


library(rgdal)
library(sp)
library(data.table)
library(ggplot2)
library(mapview)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read & prepare data -----------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read PL data (the output of the Prepare_data.R script). Select only the needed
# columns.
# ES_dt <- fread("Output/GloPL_with_id_updated_ES.csv", 
#               colClasses = "character",
#               na.strings = c("NA","N/A","null", ""),
#               select = c("unique_number", "Longitude", "Latitude"))

# Transform longitude to numeric
#ES_dt[, lon := as.numeric(Longitude)]

# Transform latitude to numeric
#ES_dt[, lat := as.numeric(Latitude)]

# Are the coordinates within expected ranges?
#range(ES_dt[,lon]) %between% c(-180, 180)
#range(ES_dt[,lat]) %between% c(-90, 90)


# Transform unprojected long-lat in Robinson coordinates
ES_dt[, c("X.prj","Y.prj") := data.table(rgdal::project(xy   = cbind(lon, lat),
                                                        proj = "+init=ESRI:54030"))]
# "+init=ESRI:54030" same as "+proj=robin"

# Check points with interactive map
points_WGS84 <- sp::SpatialPointsDataFrame(coords      = ES_dt[,.(lon, lat)], # order matters
                                           data        = ES_dt[,.(unique_number)], 
                                           proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
html_map <- mapview(points_WGS84)
# save as html
mapshot(html_map, url = "Global_map.html")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load & prepare NaturalEarth shapefiles ----------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load("Data/NaturalEarth.RData"); setDT(lbl.Y); setDT(lbl.X)
# Details about NaturalEarth shapefiles:
#   The files were already downloaded from http://www.naturalearthdata.com/
#   Graticules were adjusted to 10 dg for latitude lines and 20 dg for longitude lines 
#  (some editing was carried in ArcMap)

# Project from long-lat (unprojected) to Robinson projection
NE_countries_rob  <- spTransform(NE_countries, CRS("+proj=robin"))
NE_graticules_rob <- spTransform(NE_graticules, CRS("+proj=robin"))
NE_box_rob        <- spTransform(NE_box, CRS("+proj=robin"))

# Shift longitude of OX graticule labales. This was needed because, for example,
# the 160 dg label ended up on the 180 longitude line when projecting to
# Robinson. A shift is applied for each degree in the sequence 
# seq(from = 160, to = 0, by = -20)
shift <- c(10, 10, 9, 8, 8, 5, 2, 0, 0)
lbl.X[, shift := rep(c(shift, -rev(shift)[-1]), 2)]
lbl.X
lbl.X[, lon := lon - shift] # apply shift
lbl.X[, shift := NULL] # delete column

# Project the labales for graticules to Robinson
lbl.Y[, c("X.prj","Y.prj") := data.table(rgdal::project(xy   = cbind(lon, lat),
                                                        proj = "+proj=robin"))]
lbl.X[, c("X.prj","Y.prj") := data.table(rgdal::project(xy   = cbind(lon, lat),
                                                        proj = "+proj=robin"))]
# Create helper columns with nudged coordinates for plotting graticule labeles.
# For lbl.Y nudge longitude and for lbl.X nudge latitude.
# Give nudge values in dg (if you change, re-run also the projection lines above)
my_nudge <- cbind(nudge_lon = 10, 
                  nudge_lat = 4) 
my_nudge <- rgdal::project(my_nudge, proj = "+proj=robin")
lbl.Y[, X.prj := ifelse(lon < 0, 
                        yes = X.prj - my_nudge[1,1], 
                        no = X.prj + my_nudge[1,1])]
lbl.X[, Y.prj := ifelse(lat < 0, 
                        yes = Y.prj - my_nudge[1,2], 
                        no = Y.prj + my_nudge[1,2])]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot map ----------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ... Prepare base map ----------------------------------------------------

base_map <- 
  ggplot() +
  # ___ add graticules projected to Robinson
  geom_path(data = NE_graticules_rob, 
            aes(x     = long, 
                y     = lat, 
                group = group), 
            linetype = "dotted", 
            color    = "grey50", 
            size     = 0.1) +
  # ___ add Natural Earth countries projected to Robinson
  geom_polygon(data = NE_countries_rob, 
               aes(x     = long,
                   y     = lat, 
                   group = group), 
               colour = "grey60", # country border color
               fill   = "gray90", # country fill color
               size   = 0.2) +
  # ___ add graticule labels - latitude and longitude
  geom_text(data = lbl.Y, 
            aes(x     = X.prj, 
                y     = Y.prj, 
                label = lbl), 
            color   = "grey50", 
            size    = 1) +
  geom_text(data = lbl.X, 
            aes(x     = X.prj, 
                y     = Y.prj, 
                label = lbl), 
            color   = "grey50", 
            size    = 1) +
  # ___ add Natural Earth box projected to Robinson
  geom_polygon(data = NE_box_rob, 
               aes(x = long, 
                   y = lat), 
               colour ="black", 
               fill   ="transparent", 
               size   = 0.2) +
  # "Regions defined for each Polygons" warning has to do with fortify
  # transformation. Might get deprecated in future.
  # ___ the default ratio = 1 in coord_fixed ensures that one unit on the x-axis 
  # is the same length as one unit on the y-axis
  coord_fixed(ratio = 1) +
  # ___ remove the background and default gridlines
  theme_void()


# ... Add study locations (points) ----------------------------------------

glopl_map <- base_map +
  # ___ add the XY points
  geom_point(data = ES_dt, 
             aes(x = X.prj, 
                 y = Y.prj),
             color = "black",
             size  = 0.5,
             shape = 1,
             alpha = 1/2) +
  # Adjust theme components
  theme(
    # Set font size & family - affects legend only 
    # "sans" = "Arial" and is the default on Windows OS; check windowsFonts()
    text = element_text(size = 8, family = "sans"),
    # Grab bottom-right (x=1, y=0) legend corner 
    legend.justification = c(1, 0),
    # and position it in the bottom-right plot area.
    legend.position = c(1.05, 0.05),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    # Set height of legend items (keys).
    legend.key.height = unit(3, "mm"),
    # Set margin around entire plot.
    plot.margin = unit(c(t = 0, r = 0.8, b = 0, l = -0.5), "cm")
  )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save as pdf and png file ------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
