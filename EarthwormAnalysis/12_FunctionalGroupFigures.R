########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
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
source("Functions/FormatData.R")
source("Functions/Plots.R")
source(file.path("Functions", "ColourPicker.R"))


#################################################
# 2. Loading in variables
#################################################
models <- "Models"
data_in <-"11_Data"



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
#############################################
## models
############################################
load(file.path(models, "abundancemodel_functionalgroups.rds"))
load(file.path(models, "biomassmodel_functionalgroups.rds"))


#############################################
## BIOMASS
############################################
biomassCols <- ColourPicker(biomass$ESA)


newdata <- createNewdata(model = biomass_model, modelFixedEffects = c("ESA","variable", "scalePH", "scaleCLYPPT", "scaleSLTPPT",  
                                                                "scaleCECSOL", "bio10_12_scaled", "bio10_15_scaled", "SnowMonths_cat"),
                         mainEffect = c("ESA", "variable"), data = biomass)


## Only reference level for snow
newdata <- newdata[newdata$SnowMonths_cat == 0,]

## Predict response and variance
newdata <- predictValues(model = biomass_model, newdata, responseVar = "logValue", 
                         seMultiplier = 1.96, re.form = NA)

if(length(biomassCols) != length(levels(newdata$ESA))){
  print("The length of colours does not match the number of factor levels. Some may be removed or duplicated")
}
  
  cols <- biomassCols[1:length(levels(newdata$ESA))]
  
  cols <- rep(cols, times = length(levels(newdata$variable)))
  pt_cols <- paste("#", cols, sep="")
  
  newdata$Colour <- pt_cols
  
  #### Get rid of levels not represented by data
  fgtokeep <- c("Epi_biomass","Endo_biomass","Ane_biomass")
  newdata <- newdata[newdata$variable %in% fgtokeep,]
  
  
  ##### Group by Effect1
  order <- c()
  for(i in (levels(biomass$ESA))){
    order <- c(order, which(newdata$ESA == i))
  }
  
  labelsNew <- c("Broadleaf \n deciduous \n forest","Needleleaf \n evergreen \n forest",
                 "Mixed forest","Herbaceous \n with spare \n tree/shrub",
                 "Shrub","Herbaceous","Production - \n Herbaceous","Production - \n Plantation") 
  
  newdata <- newdata[order,]
  
  jpeg(file.path(figures, "Biomass_ESA+FG.jpg"), width = 1000, height =750, quality = 200)
  par(mar = c(5, 4, 2, 2))
  plot(-1e+05, -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
       xlim = c(0.9, nrow(newdata)),  ylab = "log Biomass", xlab = "", xaxt='n')
  
  errbar(1:nrow(newdata), newdata$logValue, newdata$upper, newdata$lower,
         add = TRUE, col = newdata$Colour, errbar.col = newdata$Colour, cex = 2,
         pch = rep(c(16, 17, 15), times = length(levels(newdata$ESA))))
  axis(1, at = c(2, 5, 8, 11, 14, 17, 20, 23), labelsNew, padj=0.5) 
  legend("topleft", legend = c("Epigeics", "Endogeics", "Anecics"), col = "Black", pch = c(16, 17, 15),  lwd = 2, bty = "n", cex = 2)
  dev.off()

  
  ##########################################
  # Abundance
  
  abundanceCols <- ColourPicker(abundance$ESA)
  
  
  newdata <- createNewdata(model = abundance_model, 
                           modelFixedEffects = c("ESA", "variable", "scalePH", "scaleCLYPPT" ,
                                                   "scaleSLTPPT",  "scaleORCDRC", "scaleCECSOL", "bio10_4_scaled", 
                                                 "bio10_15_scaled", "SnowMonths_cat","scaleAridity", "ScalePET"),
                             mainEffect = c("ESA", "variable"), data = abundance)
  
  
  ## Only reference level for snow
  newdata <- newdata[newdata$SnowMonths_cat == 0,]
  
  ## Predict response and variance
  newdata <- predictValues(model = abundance_model, newdata, responseVar = "logValue", 
                           seMultiplier = 1.96, re.form = NA)
  
  if(length(abundanceCols) != length(levels(newdata$ESA))){
    print("The length of colours does not match the number of factor levels. Some may be removed or duplicated")
  }
  
  cols <- abundanceCols[1:length(levels(newdata$ESA))]
  
  cols <- rep(cols, times = length(levels(newdata$variable)))
  pt_cols <- paste("#", cols, sep="")
  
  newdata$Colour <- pt_cols
  
  #### Get rid of levels not represented by data
  fgtokeep <- c("Epi_abundance","Endo_abundance","Ane_abundance")
  newdata <- newdata[newdata$variable %in% fgtokeep,]
  
  
  ##### Group by Effect1
  order <- c()
  for(i in (levels(abundance$ESA))){
    order <- c(order, which(newdata$ESA == i))
  }
  
  labelsNew <- c("Broadleaf \n deciduous \n forest",
                 "Broadleaf \n evergreen \n forest", "Needleleaf \n evergreen \n forest",
                 "Mixed forest","Herbaceous \n with spare \n tree/shrub",
                 "Shrub","Herbaceous","Production - \n Herbaceous","Production - \n Plantation") 
  
  newdata <- newdata[order,]
  
  jpeg(file.path(figures, "abundance_ESA+FG.jpg"), width = 1000, height =750, quality = 200)
  par(mar = c(5, 4, 2, 2))
  plot(-1e+05, -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
       xlim = c(0.9, nrow(newdata)),  ylab = "log Abundance", xlab = "", xaxt='n')
  
  errbar(1:nrow(newdata), newdata$logValue, newdata$upper, newdata$lower,
         add = TRUE, col = newdata$Colour, errbar.col = newdata$Colour, cex = 2,
         pch = rep(c(16, 17, 15), times = length(levels(newdata$ESA))))
  axis(1, at = c(2, 5, 8, 11, 14, 17, 20, 23, 26), labelsNew, padj=0.5) 
  legend("topleft", legend = c("Epigeics", "Endogeics", "Anecics"), col = "Black", pch = c(16, 17, 15),  lwd = 2, bty = "n", cex = 2)
dev.off()  



###############################################
## TRIANGULAR PLOT
##############################################
# MADE UP DATA

library(ade4)

labelsESA <- as.factor(c("Broadleaf deciduous forest",
               "Broadleaf evergreen forest", "Needleleaf evergreen forest",
               "Mixed forest","Herbaceous with spare tree/shrub",
               "Shrub","Herbaceous","Production - Herbaceous","Production - Plantation") )

Cols <- ColourPicker(labelsESA)
# This is not in the right order!!

df <- data.frame(epigeic = runif(9, 1, 3), endogeics = runif(9, 1, 3),anecics =runif(9, 1, 3))
df$total <- rowSums(df[,1:3])
df$col <- paste0("#", Cols)
row.names(df) <- labelsESA
t <- triangle.plot(df[,1:3])
points(t, col = df$col, cex = df$total, pch = 19)
# legend for colours
legend("topright", legend=labelsESA, fill = df$col, cex = 0.8, bty = "n")
# Legend for point size
pts <- seq(floor(min(df[,4])), ceiling(max(df[,4])), length.out = 3)
pts <- round(pts)
legend("bottomleft", legend= pts, fill = "black",pt.cex = pts, cex = 0.8, 
       bty = "n", pch = 19, y.intersp = 2.5, x.intersp =3)
