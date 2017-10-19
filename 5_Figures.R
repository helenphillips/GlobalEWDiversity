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

data_in <-"3_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
loadin <- files[grep(date, files)]

rm(files)
rm(date)

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

#################################################
# 3. Load in data
#################################################

sites <- read.csv(file.path(data_in, loadin))
rm(loadin)

sites$scalePH <-scale(sites$ph_new)
sites_habitat <- sites[sites$HabitatCover != "Unknown/Other",]
sites_lu <- droplevels(sites[sites$LandUse != "Unknown",])
sites <- SiteLevels(sites) 
levels(sites)[levels(sites) == "Annual crop"] <- "Annual crops"
#################################################
# 4. Load in models
#################################################

models <- "Models"


load(file.path(models, "biomass_lymgmtintensity.rds"))
load(file.path(models, "abundance_lymgmtintensity.rds"))
load(file.path(models, "abundance_a2a.rds"))


levels(b2a@frame$LU_Mgmt)[which(levels(b2a@frame$LU_Mgmt) == "Annual crop")] <- "Annual crops"
levels(a1@frame$LU_Mgmt)[which(levels(a1@frame$LU_Mgmt) == "Annual crop")] <- "Annual crops"
levels(a2a@frame$LU_Mgmt)[which(levels(a2a@frame$LU_Mgmt) == "Annual crop")] <- "Annual crops"

#################################################
# BIOMASS
#################################################
biomass <- sites[complete.cases(sites$logBiomass),]
biomass <- droplevels(biomass[biomass$LU_Mgmt != "Unknown",])
biomass$scalePH <- as.vector(scale(biomass$ph_new))




#################################################
# 5. Pick colors for figures
#################################################
# habitCols <- ColourPicker(sites_habitat$HabitatCover)
# luCols <- ColourPicker(sites_lu$LandUse)
luMgmtCols <- ColourPicker(biomass$LU_Mgmt)
#################################################
# 6. Figures
#################################################
### Interactions
#pdf(file = file.path(figures, "Biomass_ph.pdf"), height = 4)
jpeg(file = file.path(figures, "Biomass_ph.jpg"), quality = 100, res = 200, height = 1000, width = 2000)
plotSingle(model= b2a, modelFixedEffs = c("intensity", "LU_Mgmt", "scalePH"), Effect1 = "scalePH", 
                responseVar = "logBiomass", seMultiplier = 1, data = biomass, cols = "000000", 
                legend.position, ylabel = "log(Wet Biomass)", xlabel = "pH", otherContEffectsFun = "median")
dev.off()  


# plotInteraction(model = b2a, Effect1 = "scalePH", Effect2 = "LU_Mgmt", 
#                 modelFixedEffs = c("intensity", "LU_Mgmt", "scalePH"),
#                 responseVar = "logBiomass", seMultiplier = 1.96, 
#                 data = biomass, cols = luMgmtCols, legend.position = "topleft", 
#                 ylabel = "", xlabel = "")
  



# pdf(file = file.path(figures, ""), height = 4)
jpeg(file = file.path(figures, "Biomass_LuMgmt.jpg"), quality = 100, res = 200, height = 1000, width = 1500, pointsize = 12)

plotSingle(model= b2a, modelFixedEffs = c("intensity", "LU_Mgmt", "scalePH"), Effect1 = "LU_Mgmt", 
           responseVar = "logBiomass", seMultiplier = 1, data = biomass, cols = luMgmtCols, 
           legend.position, ylabel = "log(Wet Biomass)", xlabel = "", otherContEffectsFun = "median")
dev.off()


#####################################################################
## Showing intensity for one management
######################################################################

nd <- createNewdata(model= b2a, modelFixedEffects = c("intensity", "LU_Mgmt", "scalePH"), data = biomass)
nd <- predictValues(model = b2a, newdata = nd, responseVar = "logBiomass", re.form = NA, seMultiplier = 1)
nd <- nd[nd$LU_Mgmt == "Pastures (grazed lands)",]  
ref <- 0 ## Getting mean value of pH
nd <- nd[which(abs(nd$scalePH-ref)==min(abs(nd$scalePH-ref))),]

## Colours
pastureCols <- c("#cd8500","#a46a00","#7b4f00","#523500","#000000")


## Plot
# pdf(file = file.path(figures, "Biomass_PastureIntensity.pdf"), height = 4)
jpeg(file = file.path(figures, "Biomass_PastureIntensity.jpg"), quality = 100, res = 200, height = 1000, width = 1500)

par(mar=c(4.5, 4, 1, 1))
plot(-1e+05, -1e+05, ylim = c(min(nd$lower,na.rm = TRUE), max(nd$upper, na.rm = TRUE)),
     xlim = c(0, (nrow(nd)-1)),  ylab = "log(Wet Biomass)", xlab = "Management Intensity",  xaxt='n', axes = FALSE)
Axis(side = 2 )
errbar(0:(nrow(nd)-1), nd$logBiomass, nd$upper, nd$lower,
       add = TRUE, col = pastureCols, errbar.col = pastureCols, cex = 1.5)
axis(side=1, at = 0:(nrow(nd)-1), labels = levels(nd$intensity))
dev.off()


#################################################
# ABUNDANCE
#################################################
abundance <- sites[complete.cases(sites$Site_Abundance),]
abundance <- droplevels(abundance[abundance$LU_Mgmt != "Unknown",])
abundance$scalePH <- scale(abundance$ph_new)

#################################################
# 5. Pick colors for figures
#################################################
# habitCols <- ColourPicker(sites_habitat$HabitatCover)
# luCols <- ColourPicker(sites_lu$LandUse)
luMgmtCols <- ColourPicker(abundance$LU_Mgmt)
#################################################
# 6. Figures
#################################################
### Interactions
# pdf(file = file.path(figures, "Abundance_ph.pdf"), height = 4)
jpeg(file = file.path(figures, "Abundance_ph.jpg"), quality = 100, res = 200, height = 1000, width = 1500)

plotSingle(model= a2a, modelFixedEffs = c("intensity", "LU_Mgmt", "scalePH"), Effect1 = "scalePH", 
           responseVar = "logAbundance", seMultiplier = 1, data = abundance, cols = "000000", 
           legend.position, ylabel = "log(Abundance)", xlabel = "pH", otherContEffectsFun = "median")
dev.off()  




pdf(file = file.path(figures, "Abundance_LuMgmtxIntensity.pdf"), height = 4)
par(mfrow=c(1,2), mar = c(2, 2, 0, 0))
plotInteraction(model = a2a, Effect1 = "LU_Mgmt", Effect2 = "intensity", 
                modelFixedEffs = c("intensity", "LU_Mgmt", "scalePH"),
                responseVar = "logAbundance", seMultiplier = 1.96, 
                data = abundance, cols = luMgmtCols, legend.position = "topleft", 
                ylabel = "log(Abundance)", xlabel = "Management and Intensity")
cols <- paste("#", luMgmtCols, sep ="")
plot(x=(1:3), y=c(rep(2.9, 3)), xlim=c(1.9, 7), ylim = c(0,3), axes =FALSE, col = cols[1], type="l", lwd=2, ylab="", xlab="")
points(x=(1:3), y=c(rep(2.5, 3)), col = cols[2], type="l", lwd=2)
points(x=(1:3), y=c(rep(2.1, 3)), col = cols[3], type="l", lwd=2)
points(x=(1:3), y=c(rep(1.7, 3)), col = cols[4], type="l", lwd=2)
points(x=(1:3), y=c(rep(1.3, 3)), col = cols[5], type="l", lwd=2)
points(x=(1:3), y=c(rep(0.9, 3)), col = cols[6], type="l", lwd=2)
points(x=(1:3), y=c(rep(0.5, 3)), col = cols[7], type="l", lwd=2)
text(c(3.3, 3.3, 3.3), y=c(2.9,2.5,2.1,1.7,1.3,0.9,0.5), labels = c(levels(biomass$LU_Mgmt)[1:7]), adj=c(0,0.5))

dev.off()

#pdf(file = file.path(figures, "Abundance_LuMgmt.pdf"), height = 4)
jpeg(file = file.path(figures, "Abundance_LuMgmt.jpg"), quality = 100, res = 200, height = 1000, width = 1500)

plotSingle(model= a2a, modelFixedEffs = c("intensity", "LU_Mgmt", "scalePH"), Effect1 = "LU_Mgmt", 
           responseVar = "logAbundance", seMultiplier = 1, data = abundance, cols = luMgmtCols, 
           legend.position, ylabel = "log(Abundance)", xlabel = "", otherContEffectsFun = "median")
dev.off()


#####################################################################
## Showing intensity for one management
######################################################################

nd <- createNewdata(model= a2a, modelFixedEffects = c("intensity", "LU_Mgmt", "scalePH"), data = abundance)
nd <- predictValues(model = a2a, newdata = nd, responseVar = "logAbundance", re.form = NA, seMultiplier = 1)
nd <- nd[nd$LU_Mgmt == "Pastures (grazed lands)",]  
ref <- 0 ## Getting mean value of pH
nd <- nd[which(abs(nd$scalePH-ref)==min(abs(nd$scalePH-ref))),]

## Colours
pastureCols <- c("#cd8500","#a46a00","#7b4f00","#523500","#000000")


## Plot
pdf(file = file.path(figures, "Abundance_PastureIntensity.pdf"), height = 4)
par(mar=c(4.5, 4, 1, 1))
plot(-1e+05, -1e+05, ylim = c(min(nd$lower,na.rm = TRUE), max(nd$upper, na.rm = TRUE)),
     xlim = c(0, (nrow(nd)-1)),  ylab = "log(Abundance)", xlab = "Management Intensity",  xaxt='n', axes = FALSE)
Axis(side = 2 )
errbar(0:(nrow(nd)-1), nd$logAbundance, nd$upper, nd$lower,
       add = TRUE, col = "black", errbar.col = "black", cex = 1.5)
axis(side=1, at = 0:(nrow(nd)-1), labels = levels(nd$intensity))
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
