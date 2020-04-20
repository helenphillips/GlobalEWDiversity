########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\USers\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
}

#################################################
# 1. Loading libraries
#################################################
library(sp)
library(raster)
library(rgdal)

#################################################
# 2. Loading in variables
#################################################

data_in <-"2_Data"

files <- list.files(file.path(data_in))
files <- files[grep("sitesWithChelsaAndSoil_", files)]

file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinsites <- loadin[grep("sitesWithChelsaAndSoil_", loadin)]

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

sites <- read.csv(file.path(data_in, loadinsites))
# sites <- SiteLevels(sites)

#################################################
# 4. Load in TIFs
#################################################

### elevation

ele <- "I:\\sWorm\\Elevation"
tif <- raster(file.path(ele, "elevation.tif"))
sites$elevation <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees))

############## SNOW

snow <- "I:\\sWorm\\Database\\snow"
tif <- raster(file.path(snow, "snow_2015_sum.tif"))
sites$SnowMonths <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees))


############## ARIDITY
arid <- "I:\\sWorm\\Database\\AI_annual"
tif <- raster(file.path(arid, "ai_yr_TIF.tif"))
sites$Aridity <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees))
sites$scaleAridity <- scale(sites$Aridity)
############## PET

pet <- "I:\\sWorm\\Database\\PET_he_annual"
tif <- raster(file.path(pet, "pet_he_yr_TIF.tif"))
sites$PETyr <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees))
sites$ScalePET <- scale(sites$PETyr)

petSD <- "I:\\sWorm\\Database\\PET_he_monthly"
tif <- raster(file.path(petSD, "pet_he_SD.tif"))
sites$PET_SD <- extract(tif, data.frame(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees))
sites$ScalePETSD <- scale(sites$PET_SD)
rm(tif)

############# BIomes

shape <- readOGR(dsn = "I:\\sDiv\\Phillips\\sWorm\\SpatialAnalysis\\ecoregions\\official_teow\\official", layer = "wwf_terr_ecos")

test <- sites[!(is.na(sites$Latitude__decimal_degrees)),]

pts <- SpatialPoints(test[,c('Longitude__Decimal_Degrees','Latitude__decimal_degrees')], 
                     proj4string = CRS(proj4string(shape)))

test2 <- pts %over% shape

biomes <- data.frame(ID = test$ID, biome = test2$BIOME)

biomes$biome <- as.character(biomes$biome)
biomes$biome[biomes$biome ==  1] <- "Tropical & Subtropical Moist Broadleaf Forests"
biomes$biome[biomes$biome ==  2] <- "Tropical & Subtropical Dry Broadleaf Forests"
biomes$biome[biomes$biome ==  3] <- "Tropical & Subtropical Coniferous Forests"
biomes$biome[biomes$biome ==  4] <- "Temperate Broadleaf & Mixed Forests"
biomes$biome[biomes$biome ==  5] <- "Temperate Conifer Forests"
biomes$biome[biomes$biome ==  6] <- "Boreal Forests/Taig"
biomes$biome[biomes$biome ==  7] <- "Tropical & Subtropical Grasslands, Savannas & Shrublands"
biomes$biome[biomes$biome ==  8] <- "Temperate Grasslands, Savannas & Shrublands"
biomes$biome[biomes$biome ==  9] <- "Flooded Grasslands & Savannas"
biomes$biome[biomes$biome ==  10] <- "Montane Grasslands & Shrublands"
biomes$biome[biomes$biome ==  11] <- "Tundra"
biomes$biome[biomes$biome ==  12] <- "Mediterranean Forests, Woodlands & Scrub"
biomes$biome[biomes$biome ==  13] <- "Deserts & Xeric Shrublands"
biomes$biome[biomes$biome ==  14] <- "Mangroves"


sites <- merge(sites, biomes, by = 'ID', all.x = TRUE)

############## Country name

shape <- readOGR(dsn = "I:\\sDiv\\Phillips\\sWorm\\SpatialAnalysis\\world\\TM_WORLD_BORDERS_SIMPL-0.3", layer = "TM_WORLD_BORDERS_SIMPL-0.3")

test <- sites[!(is.na(sites$Latitude__decimal_degrees)),]

pts <- SpatialPoints(test[,c('Longitude__Decimal_Degrees','Latitude__decimal_degrees')], 
                     proj4string = CRS(proj4string(shape)))

test2 <- pts %over% shape
country <- data.frame(ID = test$ID, country = test2$NAME)

sites <- merge(sites, country, by = 'ID', all.x = TRUE)


###
countrytest <- sites[, c('ID', "Site_Name", 'Country', 'country')]
countrytest$Country <- as.character(countrytest$Country)
countrytest$country <- as.character(countrytest$country)

countrytest <- countrytest[which(countrytest$Country != countrytest$country),]
write.csv(countrytest, file = "C:\\Users\\hp39wasi\\temp\\checkingcountries.csv")
## I manuall checked all the sites in this file
## Any mis-matches were because of country borders (after I fixed mistakes)


#### Some sites don't have the country
nrow(sites[which(is.na(sites$country)),])
## More than just the ones missing coordinates
missingcountries <- sites[which(is.na(sites$country)),]
missingcountries <- missingcountries[,c('ID', "Site_Name", 'Country', 'country', 'Latitude__decimal_degrees', 'Longitude__Decimal_Degrees', 'bio10_1')]
write.csv(missingcountries, file = "C:\\Users\\hp39wasi\\temp\\checkingmissingcountries.csv")
## All fine, just because they are at the edge or on an island


#################################################
# 5. Save data
#################################################

write.csv(sites, file = file.path(data_out, paste("sitesWithChelsaAndSoilAndOthers_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)



