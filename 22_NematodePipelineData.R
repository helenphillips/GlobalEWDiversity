### Creating a dataframe to send to Johan van den Hooten
### Who will run it through the pipeline they developed
### for nematodes.
### Will produce another global map, but also confidence
### estimates



if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}


if(!dir.exists("DataForJohan")){
  dir.create("DataForJohan")
}
data_out <- "DataForJohan"

#################################################
# 1. Loading libraries
#################################################

source(file.path("Functions", "FormatData.R"))

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




#################################################
# 3. Load in data
#################################################


richness <- droplevels(SiteLevels(richness))
biomass <- droplevels(SiteLevels(biomass))
abundance <- droplevels(SiteLevels(abundance))

#################################################
# 4. Rename habitat cover to codes (easier for other people)
#################################################
renameESA <- function(dat){
  
  levels(dat$ESA)[levels(dat$ESA) == 'Broadleaf deciduous forest'] <- "60"
  levels(dat$ESA)[levels(dat$ESA) == 'Broadleaf evergreen forest'] <- "50"
  levels(dat$ESA)[levels(dat$ESA) == 'Needleleaf evergreen forest'] <- "70"
  levels(dat$ESA)[levels(dat$ESA) == 'Mixed forest'] <- "90"
  levels(dat$ESA)[levels(dat$ESA) == 'Herbaceous with spare tree/shrub'] <- "110"
  levels(dat$ESA)[levels(dat$ESA) == 'Shrub'] <- "120"
  levels(dat$ESA)[levels(dat$ESA) == 'Herbaceous'] <- "130"
  levels(dat$ESA)[levels(dat$ESA) == 'Production - Herbaceous'] <- "10"
  levels(dat$ESA)[levels(dat$ESA) == 'Production - Plantation'] <- "12"
  levels(dat$ESA)[levels(dat$ESA) == 'Cropland/Other vegetation mosaic'] <- "30"
  levels(dat$ESA)[levels(dat$ESA) == "Bare area (unconsolidated"] <- "202"
  
  return(dat$ESA)               
}

richness$ESA <-renameESA(richness)
abundance$ESA <-renameESA(abundance)
biomass$ESA <-renameESA(biomass)

#################################################
# 5. Format data
#################################################
neededCols <- c("ID", "Latitude__decimal_degrees", "Longitude__Decimal_Degrees", 
                "ph_new", "Clay__percent", "Silt__percent", "Organic_Carbon__percent", "ESA")


richness <- richness[,names(richness) %in% c(neededCols, "SpeciesRichness")]
richness <- richness[c(neededCols, "SpeciesRichness")]


biomass <- biomass[,names(biomass) %in% c(neededCols, "Site_Biomassm2")]
biomass <- biomass[c(neededCols, "Site_Biomassm2")]

abundance <- abundance[,names(abundance) %in% c(neededCols, "Sites_Abundancem2")]
abundance <- abundance[c(neededCols, "Sites_Abundancem2")]


#################################################
# 6. Save
#################################################


write.csv(richness, file = file.path(data_out, "sWormData_richness.csv"), row.names = FALSE)
write.csv(abundance, file = file.path(data_out, "sWormData_abundance.csv"), row.names = FALSE)
write.csv(biomass, file = file.path(data_out, "sWormData_biomass.csv"), row.names = FALSE)
