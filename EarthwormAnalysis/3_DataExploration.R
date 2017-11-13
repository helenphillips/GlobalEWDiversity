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
library(plyr)
library(dplyr)

#################################################
# 2. Loading in variables
#################################################

data_in <-"2_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

if(!dir.exists("3_Data")){
  dir.create("3_Data")
}
data_out <- "3_Data"

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"
#################################################
# 3. Load in data
#################################################

sites <- read.csv(file.path(data_in, loadin))
sites <- SiteLevels(sites) 

#################################################
# 4. Basic stats
#################################################

length(unique(sites$file)) ## 116 papers
length(unique(sites$Study_Name)) ## 143 studies

length(unique(sites$Country))## 42 Countries

#################################################
# 5. Create Map
#################################################

coord<-aggregate(cbind(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees), list(sites$Study_Name), mean)
## six don't have coordinates yet
coord <- coord[complete.cases(coord),]


coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")


# pdf(file = file.path(figures, "Map_alldata.pdf"), height = 4)
jpeg(filename = file.path(figures, "Map_alldata.jpg"), quality = 100, res = 300, width = 2000, height = 2000)
mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)
dev.off()

######################################################
## Remove Sites with no coordinates
######################################################
sites <- sites[complete.cases(sites$Latitude__decimal_degrees),]


######################################################
## 
######################################################

Summary.df <- sites %>% # Start by defining the original dataframe, AND THEN...
  group_by(Study_Name) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    N_sites = n(),  ## The function n() in dlpyr gives you the number of observations
    N_ph = sum(!(is.na(PH))),
    N_OC = sum(!(is.na(Organic_Carbon__percent))),
    N_SOM = sum(!(is.na(Soil_Organic_Matter__percent))),
    N_CN = sum(!(is.na(C.N_ratio))),
    N_Moisture = sum(!(is.na(Soil_Moisture_percent))),
    N_CEC = sum(!(is.na(CEC))),
    N_BaseS = sum(!(is.na(Base_Saturation_percent))),
    N_SandSiltClay = sum(!(is.na(Sand__percent))),
    N_STexture = sum(!(is.na(USDA_SoilTexture))),
    N_SType = sum(!(is.na(WRB_FAO_SoilType)))
  )

summary.df <- as.data.frame(Summary.df)

table(sites$LandUse)
table(sites$HabitatCover)
table(sites$HabitatCover, sites$Study_Name)
table(sites$LU_Mgmt)

###########################################################
## Looking at species richness
#############################################################


sum(!(is.na(sites$PH)))
sum(!(is.na(sites$Organic_Carbon__percent)))
sum(!(is.na(sites$Soil_Organic_Matter__percent)))
sum(!(is.na(sites$C.N_ratio)))
sum(!(is.na(sites$Soil_Moisture_percent)))
sum(!(is.na(sites$CEC)))
sum(!(is.na(sites$Base_Saturation_percent)))
sum(!(is.na(sites$Sand__percent)))
sum(!(is.na(sites$USDA_SoilTexture)))
sum(!(is.na(sites$WRB_FAO_SoilType)))

sum(is.na(sites$LandUse))
sum(sites$LandUse == "Unknown")

sum(is.na(sites$HabitatCover))
sum(sites$HabitatCover == "Unknown")

sum(is.na(sites$LU_Mgmt))
sum(sites$LU_Mgmt == "Unknown")

df.1 <- ddply(sites, .(Study_Name, Country), summarise, Variation = var(Latitude__decimal_degrees))
df.1$ Variation <- NULL
summary.df <- (merge(summary.df, df.1, by="Study_Name"))

write.csv(summary.df, file = "SoilPropertiesbyStudy.csv")


data.frame(table(sites$LandUse))
data.frame(table(sites$LU_Mgmt))
data.frame(table(sites$HabitatCover))


###########################################################
## Looking at species richness
#############################################################

hist(sites$SpeciesRichness) ## Very poisson
summary(sites$SpeciesRichness)
tapply(sites$SpeciesRichness, sites$LU_Mgmt, summary)


spR <- sites[!(is.na(sites$SpeciesRichness)),]
coord<-aggregate(cbind(spR$Longitude__Decimal_Degrees, spR$Latitude__decimal_degrees), list(spR$Study_Name), mean)
## six don't have coordinates yet
coord <- coord[complete.cases(coord),]


coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")


# pdf(file = file.path(figures, "Map_speciesRichness.pdf"), height = 4)
#jpeg(filename = file.path(figures, "Map_speciesRichness.jpg"), quality = 100, res = 300, width = 2000, height = 2000)
mar=c(0,0,0,0)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)
# dev.off()

###########################################################
## Looking at abundance
#############################################################
hist(sites$Site_Abundance)
summary(sites$Site_Abundance)
tapply(sites$Site_Abundance, sites$LU_Mgmt, summary)
sites$logAbundance <- log(sites$Site_Abundance + 1)
###########################################################
## Looking at biomass
#############################################################
hist(sites$Site_WetBiomass)
summary(sites$Site_WetBiomass)
tapply(sites$Site_WetBiomass, sites$LU_Mgmt, summary)
sites$logBiomass <- log(sites$Site_WetBiomass +1)

###########################################################
## CONTROVERSIAL!! - Altering ph values
#############################################################
sites$ph_new <- sites$PH
sites$ph_new <- ifelse(sites$PH_Collection_Method == "CaCl2", sites$PH + 1, sites$PH)
summary(sites$ph_new) # Seems reasonably after previous mistake
sites$scalePH <-scale(sites$ph_new)
###########################################################
## Distribution of the most populated soil variables
#############################################################

hist(sites$Soil_Moisture_percent)
hist(sites$Soil_Organic_Matter__percent)
hist(sites$Organic_Carbon__percent)
hist(sites$C.N_ratio)


tapply(sites$ph_new, sites$LandUse, summary)
tapply(sites$Soil_Organic_Matter__percent, sites$LandUse, summary)
tapply(sites$Organic_Carbon__percent, sites$LandUse, summary)

tapply(sites$ph_new, sites$HabitatCover, summary)
tapply(sites$Soil_Organic_Matter__percent, sites$HabitatCover, summary)
tapply(sites$Organic_Carbon__percent, sites$HabitatCover, summary)

tapply(sites$ph_new, sites$LU_Mgmt, summary)
tapply(sites$Soil_Organic_Matter__percent, sites$LU_Mgmt, summary)
tapply(sites$Organic_Carbon__percent, sites$LU_Mgmt, summary)



###########################################################
## Create a intensity variable (categorical integer)
#############################################################



intensityvars <- c("Tillage", "Pesticide", "Fertilizer",                  
                   "Selectively_harvested","Clear_cut",                   
                   "Fire","Stocking_rate",               
                   "Grazing_all_year","Rotation",                    
                   "Monoculture","Planted")


## Rotation is a good thing
sites$Rotation <- ifelse(sites$Rotation == 0, 2, sites$Rotation)
sites$Rotation <- ifelse(sites$Rotation == 1, 0, sites$Rotation)
sites$Rotation <- ifelse(sites$Rotation == 2, 1, sites$Rotation)

sites$intensity <- rowSums(sites[,which(names(sites) %in% intensityvars)], na.rm = TRUE)
sites$intensity <- ifelse(sites$LU_Mgmt == "Unknown", NA, sites$intensity)
sites$intensity <- as.factor(sites$intensity)
table(sites$LU_Mgmt, sites$intensity)
###########################################################
## Save files
#############################################################

write.csv(sites, file = file.path(data_out, paste("Sites_", Sys.Date(), ".csv", sep = "")))

######################################################
## Species level data
######################################################
###################################
speciesdata <- "0_Data"

files <- list.files(file.path(speciesdata))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadin <- loadin[grep("species_", loadin)]


#################################################
# 3. Load in data
#################################################

species <- read.csv(file.path(speciesdata, loadin))
######################################################
## How many species
######################################################
length(unique(species$SpeciesBinomial))
head(unique(species$SpeciesBinomial), 50)
