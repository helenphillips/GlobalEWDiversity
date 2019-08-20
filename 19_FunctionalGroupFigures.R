########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:/USers/hp39wasi/WORK/sWorm/EarthwormAnalysis")
}



#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
library(lme4)
library(Hmisc)
library(ade4)
source(file.path("Functions", "FormatData.R"))
source(file.path("Functions", "Plots.R"))
source(file.path("Functions", "ColourPicker.R"))
source(file.path("Functions", "cornerlabel2.R"))


plotlabcex <- 1

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
load(file.path(models, "abundancemodel_functionalgroups_revised.rds"))
load(file.path(models, "biomassmodel_functionalgroups_revised.rds"))
# load(file.path(models, "richnessmodel_functionalgroups_revised.rds"))

##########################################
# ABUNDANCE
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

jpeg(file = file.path(figures, "Biomass+Abundance_FGTriangle.jpeg"), quality = 100, res = 200, width = 2000, height = 2000)
par(mfrow=c(2, 1))
# jpeg(file = file.path(figures, "AbundanceFGTriangle.jpg"), quality = 100, res = 200, width = 2000, height = 1000)
par(mar=c(1,1,1,1))
t <- triangle.plot(df[,1:3], show.position = FALSE, min3 = c(0, 0, 0), max3 = c(1, 1, 1))
points(t, col = df$col, cex = df$total/5, pch = 19)
# legend for colours
legend(0.7,0.8, legend=labelsESA, fill = df$col, cex = 0.8, bty = "n")
corner.label2(label = "A", x = -1, y = 1, cex = plotlabcex, font = 2)

# Legend for point size
pts <- seq(floor(min(df[,5])), ceiling(max(df[,5])), length.out = 4)
pts <- round(pts)
legend(-1.5,0.8, legend= pts, pt.cex = pts/5, cex = 0.8, 
       bty = "n", pch = rep(19, times = 4), y.intersp = 2.7, x.intersp =3)
mtext(expression("Total Abundance (ind. per" ~ m^{2} ~ ")"), 3, -3, adj=0.07)

# dev.off()


#############################################
## BIOMASS
############################################
biomassCols <- ColourPicker(biomass$ESA)


newdata <- createNewdata(model = biomass_model, modelFixedEffects = c("ESA","variable", "ScaleElevation",  "scalePH", "scaleCLYPPT", "scaleSLTPPT",  "scaleORCDRC",
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



# jpeg(file = file.path(figures, "BiomassFGTriangle.jpeg"), quality = 100, res = 200, width = 2000, height = 1000)
par(mar=c(1, 1, 1, 1))
t <- triangle.plot(df[,1:3], show.position = FALSE, min3 = c(0, 0, 0), max3 = c(1, 1, 1), cpoint = 0)
points(t, col = df$col, cex = df$total/5 + 1, pch = 19)
# Legend for point size
pts <- seq(ceiling(min(df[,5])), ceiling(max(df[,5])), length.out = 4)
pts <- round(pts)
corner.label2(label = "B", x = -1, y = 1, cex = plotlabcex, font = 2)

legend(-1.5,0.8, legend= pts, pt.cex = pts/5 + 1, cex = 0.8, 
       bty = "n", pch = rep(19, times = 4), y.intersp = 2.0, x.intersp =3)
mtext(expression("Total Biomass (grams per" ~ m^{2} ~ ")"), 3, -3, adj=0.08)


# legend for colours
legend(0.7,0.8, legend=labelsESA, fill = df$col, cex = 0.8, bty = "n")
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

#################################
################# Mean plots
############################

jpeg(file = file.path(figures, "HabitatCover_FGAbundance+Biomass.jpg"), quality = 100, res = 200, width = 1800, height = 2000)
par(mfrow = c(2, 1))

##############################
## ABUNDANCE
##############################
# jpeg(file = file.path(figures, "HabitatCover_Richness+Abundance.jpg"), quality = 100, res = 200, width = 1500, height = 2000)

abundanceCols <- ColourPicker(abundance$ESA)

newdata <- createNewdata(model = abundance_model, modelFixedEffects = c("ESA","variable", "scalePH", "scaleCLYPPT", "scaleSLTPPT",  
                                                                      "scaleORCDRC", "scaleCECSOL", "bio10_7_scaled" , "scaleAridity",
                                                                      "bio10_15_scaled", "SnowMonths_cat", "ScalePET"),
                         mainEffect = c("ESA", "variable"), data = biomass)

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


newdata <- newdata[order,]
# jpeg(file.path(figures, "Biomass_ESA+FG.jpg"), width = 1000, height =750, quality = 200)
par(mar = c(6, 4, 2, 2))
plot(-1e+05, -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
     xlim = c(0.9, nrow(newdata)),  ylab = "(ln-) Abundance", xlab = "", xaxt='n', axes = FALSE)
Axis(side = 2, cex.axis = 1)

xpos <- 1:nrow(newdata)
adj1 <- c(0.2, 0, -0.2)
xpos <- xpos + adj1

errbar(xpos, newdata$logValue, newdata$upper, newdata$lower,
       add = TRUE, col = newdata$Colour, errbar.col = newdata$Colour, cex = 2,
       pch = rep(c(16, 17, 15), times = length(levels(newdata$ESA))))

mtext("(a)", side = 3, line = 0, at = 0, adj = 0.1)

labelsNew <- c("Broadleaf\ndeciduous\nforest", "Broadleaf\nevergreen\nforest", "Needleleaf\nevergreen\nforest",
               "Mixed forest","Herbaceous\nwith spare\ntree/shrub",
               "Shrub","Herbaceous","Production -\nHerbaceous","Production -\nPlantation") 
axis(1, at = c(2, 5, 8, 11, 14, 17, 20, 23, 26), labelsNew, padj=0.5, las = 2) 


##############################
## BIOMASS
##############################



biomassCols <- ColourPicker(biomass$ESA)

newdata <- createNewdata(model = biomass_model, modelFixedEffects = c("ESA","variable", "ScaleElevation",  "scalePH", "scaleCLYPPT", "scaleSLTPPT",  
                                                                      "scaleORCDRC", "scaleCECSOL", "bio10_4_scaled" , "bio10_12_scaled", 
                                                                      "bio10_15_scaled", "SnowMonths_cat", "ScalePETSD"),
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


newdata <- newdata[order,]
# jpeg(file.path(figures, "Biomass_ESA+FG.jpg"), width = 1000, height =750, quality = 200)
par(mar = c(6, 4, 2, 2))
plot(-1e+05, -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
     xlim = c(0.9, nrow(newdata)),  ylab = "(ln-) Biomass", xlab = "", xaxt='n',axes = FALSE)
Axis(side = 2, cex.axis = 1)

mtext("(b)", side = 3, line = 0, at = 0, adj = 0.1)


xpos <- 1:nrow(newdata)
adj1 <- c(0.2, 0, -0.2)
xpos <- xpos + adj1

errbar(xpos, newdata$logValue, newdata$upper, newdata$lower,
       add = TRUE, col = newdata$Colour, errbar.col = newdata$Colour, cex = 2,
       pch = rep(c(16, 17, 15), times = length(levels(newdata$ESA))))


labelsNew <- c("Broadleaf\ndeciduous\nforest", "Needleleaf\nevergreen\nforest",
               "Mixed forest","Herbaceous\nwith spare\ntree/shrub",
               "Shrub","Herbaceous","Production -\nHerbaceous","Production -\nPlantation") 
axis(1, at = c(2, 5, 8, 11, 14, 17, 20, 23), labelsNew, padj=0.5, las = 2) 



# legend("topleft", legend = c("Epigeics", "Endogeics", "Anecics"), col = "Black", pch = c(16, 17, 15),  lwd = 2, bty = "n", cex = 2)
# dev.off()





# legend("topleft", legend = c("Epigeics", "Endogeics", "Anecics"), col = "Black", pch = c(16, 17, 15),  lwd = 2, bty = "n", cex = 2)
 dev.off()

