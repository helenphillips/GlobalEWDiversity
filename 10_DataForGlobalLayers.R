## SET WD ----------------------------
setwd("C:/restore2/hp39wasi/sWorm/EarthwormAnalysis")

## FUNCTIONS AND LIBRARIES ----------------

source("Functions/loadMostRecent.R")

data_out <- "10_Data"

if(!dir.exists(data_out)){
  dir.create(data_out)
}

## LOAD THE DATA -------------------------

data_in <- "8_Data"

print("Loading datasets")
richness <- loadMostRecent("sitesRichness", data_in)
richness <- read.csv(file.path(data_in, richness ))
abundance <- loadMostRecent("sitesAbundance_", data_in)
abundance <- read.csv(file.path(data_in, abundance))
biomass <- loadMostRecent("sitesBiomass_", data_in)
biomass <- read.csv(file.path(data_in, biomass))

## ONLY KEEP RELEVANT COLUMNS --------------------

neededCols <- c("ID", "Latitude__decimal_degrees", "Longitude__Decimal_Degrees", 
                "ESA", "SnowMonths_cat",              
                "phFinal","ClayFinal","SiltFinal", "OCFinal", "CECSOL", 
                "bio10_1", "bio10_4","bio10_7",
                "bio10_12","bio10_15",
                "PETyr", "PET_SD", "Aridity", 
                "elevation" )

richness <- richness[names(richness) %in% neededCols]
abundance <- abundance[names(abundance) %in% neededCols]
biomass <- biomass[names(biomass) %in% neededCols]

## RENAMING ESA LAYERS ---------------------------------

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


richness$ESA <- renameESA(richness)
abundance$ESA <- renameESA(abundance)
biomass$ESA <- renameESA(biomass)

## SORTING SNOW VALUES ----------------------------

levels(richness$SnowMonths_cat)[which(levels(richness$SnowMonths_cat) == "4plus")] <- "4"
levels(abundance$SnowMonths_cat)[which(levels(abundance$SnowMonths_cat) == "4plus")] <- "4"
levels(biomass$SnowMonths_cat)[which(levels(biomass$SnowMonths_cat) == "4plus")] <- "4"


## RENAME TO MATCH LAYERS ----------------------------

names(richness) == names(biomass)
names(biomass) == names(abundance)

newNames <- c("ID", "Latitude__decimal_degrees", "Longitude__Decimal_Degrees", 
              "CHELSA_bio10_1", "CHELSA_bio10_4","CHELSA_bio10_7",
              "CHELSA_bio10_12","CHELSA_bio10_15",
              "CECSOL_weighted", "elevation",
              "ai_yr_TIF", "pet_he_yr_TIF", "pet_he_SD",  
              "ESACCI-LC-L4-LCCS-Map-1000m-P5Y-2010-v1.6.1_NEW", "snow_2015_sum",         
              "PHIHOX_weighted","CLYPPT_weighted","SLTPPT_weighted", "ORCDRC_weighted")


names(richness) <- names(biomass) <- names(abundance) <- newNames


write.csv(richness, file = file.path(data_out, paste0("sitesRichness_", Sys.Date(), ".csv")), row.names = FALSE)
write.csv(abundance, file = file.path(data_out, paste0("sitesAbundance_", Sys.Date(), ".csv")), row.names = FALSE)
write.csv(biomass, file = file.path(data_out, paste0("sitesBiomass_", Sys.Date(), ".csv")), row.names = FALSE)
