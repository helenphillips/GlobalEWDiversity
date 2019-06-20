########################################################
# 0. Set Working Directory
########################################################
 
if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
library(plyr)
library(dplyr)
source("Functions/FormatData.R")
#################################################
# 2. Loading in variables
#################################################

data_in <-"5_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

if(!dir.exists("6_Data")){
  dir.create("6_Data")
}
data_out <- "6_Data"

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
# 3.5 Old data  
#################################################
# During the revision stage we re-ran analysis
# But as data had still been added to teh database
# we don't want to add any new data
# but things were corrected (coordinates) so we still
# took a new cut of the database, then subset
# by the studies in the old

old_sites <- read.csv(file.path(data_in, "sites_2018-09-24.csv"))
old_studies <- unique(old_sites$file)


sites <- sites[sites$file %in% old_studies,]
#################################################
# 4. Basic stats
#################################################

length(unique(sites$file)) ## 196 papers
length(unique(sites$Study_Name)) ##  250 studies

length(unique(sites$Country))## 65 Countries
nrow(sites) # 7805
#################################################
# 5. Create Map
#################################################

coord<-aggregate(cbind(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees), list(sites$Study_Name), mean)
## 3 don't have coordinates yet
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
## Coordinates sanity check
######################################################

any(sites$Latitude__decimal_degrees == sites$Longitude__Decimal_Degrees) ## Apparently a common mistake
plot(sites$Sites_Abundancem2 ~ sites$bio10_1)
plot(sites$Sites_Abundancem2 ~ sites$bio10_7)
plot(sites$Sites_Abundancem2 ~ sites$SnowMonths)

plot(sites$Site_Biomassm2 ~ sites$bio10_1)
plot(sites$Site_Biomassm2 ~ sites$bio10_7)
plot(sites$Site_Biomassm2 ~ sites$SnowMonths)

## But there is one study where we identified that some coordiantes may be wrong
## Meta-data says they are from germany
## But original coordinates are different

sites <- sites[-(which(sites$Study_Name == "birkhofer2012" & sites$country != "Germany")),] #  7699
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
hist(sites$Sites_Abundancem2)
summary(sites$Sites_Abundancem2)
tapply(sites$Sites_Abundancem2, sites$LU_Mgmt, summary)
sites$logAbundance <- log(sites$Sites_Abundancem2 + 1)
hist(sites$logAbundance)
###########################################################
## Looking at biomass
#############################################################
hist(sites$Site_Biomassm2)
summary(sites$Site_Biomassm2)
tapply(sites$Site_Biomassm2, sites$LU_Mgmt, summary)
sites$logBiomass <- log(sites$Site_Biomassm2 +1)
hist(sites$logBiomass)
###########################################################
## CONTROVERSIAL!! - Altering ph values
#############################################################
table(sites$PH_Collection_Method)


sites$ph_new <- sites$PH
sites$ph_new[which(sites$PH_Collection_Method == "CaCl2")] <- sites$PH[which(sites$PH_Collection_Method == "CaCl2")] - 1
sites$ph_new[which(sites$PH_Collection_Method == "KCl")] <- sites$PH[which(sites$PH_Collection_Method == "KCl")] + 1
summary(sites$ph_new) # Seems reasonably after previous mistake
sites$scalePH <-scale(sites$ph_new)
###########################################################
## Distribution of the most populated soil variables
#############################################################

hist(sites$Soil_Moisture_percent)
## Not sure its the right measurement
sites$Soil_Moisture_percent[which(sites$file == "4714_Birkhofer2011")] <- NA
hist(sites$Soil_Moisture_percent)

hist(sites$Soil_Organic_Matter__percent)
## According to 4336_Scharenbroch2012 there are values above 100%

hist(sites$Organic_Carbon__percent)
## One study is obviously not the right units
sites$Organic_Carbon__percent[which(sites$file == "4327_Wu2012")] <- NA
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


##########################################################
## Snow months to categorical variable
#############################################################

## Grouping anything 5 and above into one category will give around 272 sites

sites$SnowMonths_cat <- (sites$SnowMonths)
sites$SnowMonths_cat[which(sites$SnowMonths_cat > 3)] <- 4
sites$SnowMonths_cat <- as.factor(sites$SnowMonths_cat)
levels(sites$SnowMonths_cat)[levels(sites$SnowMonths_cat) == 4] <- "4plus"


###########################################################
## Save files
#############################################################

write.csv(sites, file = file.path(data_out, paste("Sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
