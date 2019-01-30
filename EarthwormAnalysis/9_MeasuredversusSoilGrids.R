########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
library(lme4)
library(plyr)
library(dplyr)
source("Functions/FormatData.R")
#################################################
# 2. Loading in variables
#################################################

data_in <-"3_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

if(!dir.exists("3.5_Data")){
  dir.create("3.5_Data")
}
data_out <- "3.5_Data"

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
# PH
#################################################
# scalePH  + 

ph1 <- lm(ph_new ~ PHIHOX, data = sites)
plot(ph1)
plot(ph_new ~ PHIHOX, data = sites)
# r2 = 0.23

## Just sites where coordinates varied within a study 
## (we know that if coordiantes don't vary the fit won't be good!)


Summary.df <- sites %>% # Start by defining the original dataframe, AND THEN...
  group_by(Study_Name) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    coords = var(Longitude__Decimal_Degrees),
    ph = var(ph_new),
    OC = var(Organic_Carbon__percent),
    SOM = var(Soil_Organic_Matter__percent),
    CN = var(C.N_ratio),
    Moisture = var(Soil_Moisture_percent),
    Sand = var(Sand__percent),
    Silt = var(Silt__percent),
    Clay = var(Clay__percent)
  )

summary.df <- as.data.frame(Summary.df)
varStudyes <- summary.df$Study_Name[which(summary.df$coords > 0)] 


vardat <- sites[sites$Study_Name %in% varStudyes,]
ph2 <- lm(ph_new ~ PHIHOX, data = vardat)
plot(ph2)
plot(ph_new ~ PHIHOX, data = vardat)
## r2 = 0.35


###### What about the mean of each study
mean.df <- sites %>% # Start by defining the original dataframe, AND THEN...
  group_by(Study_Name) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    # coords = var(Longitude__Decimal_Degrees),
    ph = mean(ph_new),
    sgph = mean(PHIHOX)
   # OC = var(Organic_Carbon__percent),
    # SOM = var(Soil_Organic_Matter__percent),
    # CN = var(C.N_ratio),
    # Moisture = var(Soil_Moisture_percent),
    # SandSiltClay = var(Sand__percent)
  )

mean.df <- as.data.frame(mean.df)

ph3 <- lm(ph ~ sgph, data = mean.df)
plot(ph3)
plot(ph ~ sgph, data = mean.df)
# r2 = 0.392

### Regardless, stick with the approach of filling in

## 
sites$phFinal <- sites$ph_new
## Where there are NAs, use SoilGrids
sites$phFinal[which(is.na(sites$phFinal))] <- sites$PHIHOX[which(is.na(sites$phFinal))]
## Where coordinates vary but ph stays the same
t <- summary.df$Study_Name[intersect(which(summary.df$coords > 0), which(summary.df$ph == 0))]
sites$phFinal[which(sites$Study_Name %in% t)] <- sites$PHIHOX[which(sites$Study_Name %in% t)]

summary(sites$phFinal)

#################################################
# CLAY/SAND/SILT
#################################################
# scaleCLYPPT + scaleSLTPPT + 

sites$ClayFinal <- sites$Clay__percent
## Where there are NAs, use SoilGrids
sites$ClayFinal[which(is.na(sites$ClayFinal))] <- sites$CLYPPT[which(is.na(sites$ClayFinal))]
## Where coordinates vary but ph stays the same
t <- summary.df$Study_Name[intersect(which(summary.df$coords > 0), which(summary.df$Clay == 0))]
# none match criteria
# sites$ClayFinal[which(sites$Study_Name %in% t)] <- sites$CLYPPT[which(sites$Study_Name %in% t)]


sites$SandFinal <- sites$Sand__percent
## Where there are NAs, use SoilGrids
sites$SandFinal[which(is.na(sites$SandFinal))] <- sites$SNDPPT[which(is.na(sites$SandFinal))]
## Where coordinates vary but ph stays the same
t <- summary.df$Study_Name[intersect(which(summary.df$coords > 0), which(summary.df$Sand == 0))]
# none match criteria
# sites$SandFinal[which(sites$Study_Name %in% t)] <- sites$SNDPPT[which(sites$Study_Name %in% t)]

sites$SiltFinal <- sites$Silt__percent
## Where there are NAs, use SoilGrids
sites$SiltFinal[which(is.na(sites$SiltFinal))] <- sites$SLTPPT[which(is.na(sites$SiltFinal))]
## Where coordinates vary but ph stays the same
t <- summary.df$Study_Name[intersect(which(summary.df$coords > 0), which(summary.df$Silt == 0))]
# none match criteria
# sites$SiltFinal[which(sites$Study_Name %in% t)] <- sites$SLTPPT[which(sites$Study_Name %in% t)]

#################################################
# CEC
#################################################
# scaleCECSOL 

## Not worth doing with only a few hundred values in many different units

#################################################
# ORGANIC CARBON
#################################################
# scaleORCDRC

sites$OCFinal <- sites$Organic_Carbon__percent
## Where there are NAs, use SoilGrids
sites$OCFinal[which(is.na(sites$OCFinal))] <- sites$ORCDRC[which(is.na(sites$OCFinal))]
## Where coordinates vary but ph stays the same
t <- summary.df$Study_Name[intersect(which(summary.df$coords > 0), which(summary.df$OC == 0))]
# none match criteria
# sites$OCFinal[which(sites$Study_Name %in% t)] <- sites$ORCDRC[which(sites$Study_Name %in% t)]


###########################################################
## Save files
#############################################################

write.csv(sites, file = file.path(data_out, paste("Sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

