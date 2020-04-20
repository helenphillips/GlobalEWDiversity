########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}


if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\USers\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
}
source("Functions/cornerlabel2.R")

source("Functions/FormatData.R")
source("Functions/createMap.R")
library(maps)
library(maptools)

########################################################
# 2. Create folder if it doesn't exist to save data into
########################################################

if(!dir.exists("5_Data")){
  dir.create("5_Data")
}
data_out <- "5_Data"

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

#################################################
# 4 Load in new data (contains all zeros)
#################################################

# just the data that was used in one of the three corrected models
data_folder <- "8_Data"

new_richness <- read.csv(file.path(data_folder, "sitesRichness_2019-11-28.csv"))
new_abundance <- read.csv(file.path(data_folder, "sitesAbundance_2019-11-25.csv"))
new_biomass <- read.csv(file.path(data_folder, "sitesBiomass_2019-11-25.csv"))


#################################################
# 4.5 Load in older/missing-zero data
#################################################

# just the data that was used in (one) of the three original models

richness <- read.csv(file.path(data_folder, "sitesRichness_2019-06-20.csv"))
abundance <- read.csv(file.path(data_folder, "sitesAbundance_2019-06-20.csv"))
biomass <- read.csv(file.path(data_folder, "sitesBiomass_2019-06-20.csv"))


richness$studyID <- paste(richness$file, richness$Study_Name, sep = "_")
abundance$studyID <- paste(abundance$file, abundance$Study_Name, sep = "_")
biomass$studyID <- paste(biomass$file, biomass$Study_Name, sep = "_")

#################################################3
## COMPARISON MAP 
###################################################3
## Doubling of sites - Richness
oldD <- data.frame(table(richness$studyID))
newD <- data.frame(table(new_richness$studyID))
compareRichness <- merge(oldD, newD, by = "Var1", all = TRUE)
names(compareRichness) <- c("studyID", "oldN", "newN")
compareRichness$diff <- compareRichness$newN - compareRichness$oldN
compareRichness$percentChange <- (compareRichness$diff / compareRichness$oldN) * 100

doubled <- compareRichness[compareRichness$percentChange > 100,]
doubled <- doubled[complete.cases(doubled),]

## How many studies were affected
n_affected <-compareRichness[compareRichness$diff > 0,]
n_affected <- n_affected[complete.cases(n_affected),]
nrow(n_affected) # 64 studies - this doesn't include the studies that have been removed following even more data cleaning

new_richness <- droplevels(new_richness)
oldD <- data.frame(table(richness$file))
newD <- data.frame(table(new_richness$file))
compareRichness_file <- merge(oldD, newD, by = "Var1", all = TRUE)
names(compareRichness_file) <- c("studyID", "oldN", "newN")
compareRichness_file$diff <- compareRichness_file$newN - compareRichness_file$oldN

## How many files were affected
n_affected <-compareRichness_file[compareRichness_file$diff > 0,]
n_affected <- n_affected[complete.cases(n_affected),]

nrow(n_affected) # 51 files




## Doubling of sites - Abundance
oldD <- data.frame(table(abundance$studyID))
newD <- data.frame(table(new_abundance$studyID))
compareAbundance <- merge(oldD, newD, by = "Var1", all = TRUE)
names(compareAbundance) <- c("studyID", "oldN", "newN")
compareAbundance$diff <- compareAbundance$newN - compareAbundance$oldN
compareAbundance$percentChange <- (compareAbundance$diff / compareAbundance$oldN) * 100

doubled <- compareAbundance[compareAbundance$percentChange > 100,]
doubled <- doubled[complete.cases(doubled),]
nrow(doubled) # 14 doubled in size

## How many studies were affected
n_affected <-compareAbundance[compareAbundance$diff > 0,]
n_affected <- n_affected[complete.cases(n_affected),]

nrow(n_affected) # 69 studies

new_abundance <- droplevels(new_abundance)
oldD <- data.frame(table(abundance$file))
newD <- data.frame(table(new_abundance$file))
compareAbundance_file <- merge(oldD, newD, by = "Var1", all = TRUE)
names(compareAbundance_file) <- c("studyID", "oldN", "newN")
compareAbundance_file$diff <- compareAbundance_file$newN - compareAbundance_file$oldN

## How many files were affected
n_affected <-compareAbundance_file[compareAbundance_file$diff > 0,]
n_affected <- n_affected[complete.cases(n_affected),]

nrow(n_affected) # 56 files




## Doubling of sites - Biomass
oldD <- data.frame(table(biomass$studyID))
newD <- data.frame(table(new_biomass$studyID))
compareBiomass <- merge(oldD, newD, by = "Var1", all = TRUE)
names(compareBiomass) <- c("studyID", "oldN", "newN")
compareBiomass$diff <- compareBiomass$newN - compareBiomass$oldN
compareBiomass$percentChange <- (compareBiomass$diff / compareBiomass$oldN) * 100

doubled <- compareBiomass[compareBiomass$percentChange > 100,]
## There's a couple of rows, that are NAs because the studies were added - they seem fine.
## When the study name changed, for those studies, there were no change in teh number of sites
## So still ok to remove the NAs
doubled <- doubled[complete.cases(doubled),]
nrow(doubled) # 7
## How many studies were affected
n_affected <-compareBiomass[compareBiomass$diff > 0,]
n_affected <- n_affected[complete.cases(n_affected),]

nrow(n_affected) # 36 studies


new_biomass <- droplevels(new_biomass)
oldD <- data.frame(table(biomass$file))
newD <- data.frame(table(new_biomass$file))
compareBiomass_file <- merge(oldD, newD, by = "Var1", all = TRUE)
names(compareBiomass_file) <- c("studyID", "oldN", "newN")
compareBiomass_file$diff <- compareBiomass_file$newN - compareBiomass_file$oldN

## How many files were affected
n_affected <-compareBiomass_file[compareBiomass_file$diff > 0,]
n_affected <- n_affected[complete.cases(n_affected),]

nrow(n_affected) #  30 files



## In total how many files and studies
length(unique(new_richness$file)) # 123
length(unique(new_richness$studyID)) # 161

length(unique(new_abundance$file)) # 162
length(unique(new_abundance$studyID)) # 207
length(unique(new_biomass$file)) # 98
length(unique(new_biomass$studyID)) # 128



names(new_richness)[names(new_richness) == "scaleElevation"] <- "ScaleElevation"

new_richness <- droplevels(new_richness)
new_abundance <- droplevels(new_abundance)
new_biomass <- droplevels(new_biomass)

alldat <- rbind(new_richness, new_abundance, new_biomass)
length(unique(alldat$file)) # 176
length(unique(alldat$Study_Name)) # 224



## Change the name of one study which used to be NA
abundance$file[is.na(abundance$Study_Name)]
abundance$Study_Name <- as.character(abundance$Study_Name)
abundance$Study_Name[is.na(abundance$Study_Name)] <- "hurisso"
abundance$Study_Name <- as.factor(abundance$Study_Name) # Urgh hate myself for that
abundance$Study_Name[grep("Hurisso2011", abundance$file)] # Sanity check

abundance$Study_site <- as.character(abundance$Study_site)

abundance$Study_site[abundance$file == "4836_Hurisso2011"] <- 
  gsub("NA", "hurisso", abundance$Study_site[abundance$file == "4836_Hurisso2011"])

abundance$Study_site <- as.factor(abundance$Study_site) # Urgh hate myself for that

new_abundance$Study_site[grep("Hurisso2011", new_abundance$file)]




biomass$file[is.na(biomass$Study_Name)]
biomass$Study_Name <- as.character(biomass$Study_Name)
biomass$Study_Name[is.na(biomass$Study_Name)] <- "hurisso"
biomass$Study_Name <- as.factor(biomass$Study_Name) # Urgh hate myself for that
biomass$Study_Name[grep("Hurisso2011", biomass$file)] # Sanity check

biomass$Study_site <- as.character(biomass$Study_site)

biomass$Study_site[biomass$file == "4836_Hurisso2011"] <- 
  gsub("NA", "hurisso", biomass$Study_site[biomass$file == "4836_Hurisso2011"])

biomass$Study_site <- as.factor(biomass$Study_site) # Urgh hate myself for that

new_biomass$Study_site[grep("Hurisso2011", new_biomass$file)]



## Remove sites that were previously used - leaving just the zeros 
new_richness <- droplevels(new_richness[!(new_richness$Study_site %in% unique(richness$Study_site)),])
new_abundance <- droplevels(new_abundance[!(new_abundance$Study_site %in% unique(abundance$Study_site)),])
new_biomass <- droplevels(new_biomass[!(new_biomass$Study_site %in% unique(biomass$Study_site)),])

length(unique(new_richness$file)) # 52
length(unique(new_richness$studyID)) # 65
# Issue is one file/study had a decrease in the number of sites, unrelated to the error



length(unique(new_abundance$file)) # 57
length(unique(new_abundance$studyID)) # 70
 
length(unique(new_biomass$file)) # 31
length(unique(new_biomass$studyID)) # 39

## These numbers don't match previous calculations
# Manualy checking
# they are all because of changes in names or new studies 
# (unrelated to the error, but should have always  been there)



new_richness <- droplevels(new_richness[-(grep("Kernecker", new_richness$file)),])

new_abundance <- droplevels(new_abundance[-(grep("Kernecker", new_abundance$file)),])

new_biomass <- droplevels(new_biomass[-(grep("113_Davalos2015_regional", new_biomass$studyID)),])
new_biomass <- droplevels(new_biomass[-(grep("3331_Watmough2014_Watmouth2014", new_biomass$studyID)),])
new_biomass <- droplevels(new_biomass[-(grep("5597_Hirth2009_hirth1994", new_biomass$studyID)),])


## remove NAs from the two studies that have a couple of coordinates missing
new_richness <- new_richness[!(is.na(new_richness$Latitude__decimal_degrees)),]
new_abundance <- new_abundance[!(is.na(new_abundance$Latitude__decimal_degrees)),]
new_biomass <- new_biomass[!(is.na(new_biomass$Latitude__decimal_degrees)),]


wide_inch <- 4.75
point_size <- 7
plotlabcex <- 1
#legendcex <- 0.9
# resdpi <- 300

#png(file = file.path(figures, "Maps+newZeroData_correction.png"), res = 300, width = 1000, height = 3000)
# png(file.path(figures, "Maps+newZeroData_correction.png"),width=(3 * 17.5),height=(3*8.75),units="cm",res=300)
pdf(file.path(figures, "Maps+newZeroData_correction.pdf"),width= wide_inch + 2, height= (wide_inch), pointsize = point_size)

par(oma = c(0, 0, 0, 0))

par(mar = c(1, 1, 1, 1))
par(mfrow = c(3,2))
createSizedMap(richness)
corner.label2(label = "A", x = -1, y = 1, cex=plotlabcex, font = 2)
createSizedMap(new_richness)
corner.label2(label = "B", x = -1, y = 1, cex=plotlabcex, font = 2)
createSizedMap(abundance)
corner.label2(label = "C", x = -1, y = 1, cex=plotlabcex, font = 2)
createSizedMap(new_abundance)
corner.label2(label = "D", x = -1, y = 1, cex=plotlabcex, font = 2)
createSizedMap(biomass)
corner.label2(label = "E", x = -1, y = 1, cex=plotlabcex, font = 2)
createSizedMap(new_biomass)
corner.label2(label = "F", x = -1, y = 1, cex=plotlabcex, font = 2)

dev.off()
