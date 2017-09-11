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
library(dplyr)
#################################################
# 2. Loading in variables
#################################################

data_in <-"2_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
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

#################################################
# 4. Basic stats
#################################################

length(unique(sites$file)) ## 61 papers
length(unique(sites$Study_Name)) ## 78 studies

length(unique(sites$Country))## 31 Countries

#################################################
# 5. Create Map
#################################################

coord<-aggregate(cbind(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees), list(sites$Study_Name), mean)
## Five don't have coordinates yet
coord <- coord[complete.cases(coord),]


coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")


pdf(file = file.path(figures, "Map.pdf"), height = 4)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)
dev.off()



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
    N_Moisture = sum(!(is.na(Soil_Moisture_percent)))
  )

summary.df <- as.data.frame(Summary.df)

table(sites$LandUse)
table(sites$HabitatCover)
table(sites$HabitatCover, sites$Study_Name)
table(sites$LU_Mgmt)

###########################################################
## Looking at species richness
#############################################################

hist(sites$SpeciesRichness) ## Very poisson
summary(sites$SpeciesRichness)

###########################################################
## Looking at abundance
#############################################################
hist(sites$Site_Abundance)
summary(sites$Site_Abundance)

###########################################################
## Looking at species richness
#############################################################
hist(sites$Site_WetBiomass)
summary(sites$Site_WetBiomass)


###########################################################
## CONTROVERSIAL!! - Altering ph values
#############################################################
sites$ph_new <- sites$PH
sites$ph_new <- ifelse(sites$PH_Collection_Method == "CaCl2", sites$PH + 1, sites$PH)
summary(sites$ph_new) # Seems reasonably after previous mistake

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

write.csv(sites, file = file.path(data_out, paste("Sites_", Sys.Date(), ".csv", sep = "")))
