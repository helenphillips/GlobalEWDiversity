
########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

source("Functions/FormatData.R")

library(dplyr)
########################################################
# 2. Create folder if it doesn't exist to save data into
########################################################

if(!dir.exists("10_Data")){
  dir.create("10_Data")
}
data_out <- "10_Data"


#################################################
# 3. Loading in variables
#################################################
## Species data

data_in_spp <-"9_Data"
files <- list.files(file.path(data_in_spp))
files <- files[grep("UniqueSpeciestoSend_", files)]
## Change this when I get a file back from the EW experts

file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinfg <- loadin[grep("UniqueSpeciestoSend_", loadin)]

## Site data


data_in_sites <-"3.5_Data"

files <- list.files(file.path(data_in_sites))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

rm(files)
rm(date)

### Species data

data_in_spp2 <-"0_Data"
files <- list.files(file.path(data_in_spp2))
files <- files[grep("species_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin_spp <- files[grep(date, files)]

#################################################
# 4. Load in data
#################################################

sites <- read.csv(file.path(data_in_sites, loadin))
spp <- read.csv(file.path(data_in_spp, loadinfg))
spp_dat <- read.csv(file.path(data_in_spp2, loadin_spp))

#################################################
# 5. Tidy up categories in spp
#################################################

levels(spp$Revised_fg)
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "anecic"] <- "Anecic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "endogeic"] <- "Endogeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "endogeic?"] <- "Endogeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "Epi-endogeic"] <- "Epi-Endogeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "epigeic"] <- "Epigeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "N/A"] <- "Unknown"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "endogeic or epiendogeic"] <- "Epi-Endogeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "Endo-anecic"] <- "Anecic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "Endogeic, Anecic"] <- "Anecic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "epi-anecic"] <- "Anecic"

#################################################
# 6. Remove unwanted columns
#################################################
keep <- c("SpeciesBinomial","Revised", "Revised_fg")      

spp <- spp[,names(spp) %in% keep]

#################################################
# 7. Merge with species level dataset
#################################################

spp_dat <- merge(spp_dat, spp, by.x = "SpeciesBinomial", by.y = "SpeciesBinomial", all.x = TRUE)

#################################################
# 8. Use FG given if no other (especially as some non-species will have them)
#################################################

missing <- which(is.na(spp_dat$Revised_fg))
given <- which(!(is.na(spp_dat$Functional_Type)))
toFill <- intersect(missing, given)

spp_dat$Revised_fg[toFill] <- spp_dat$Functional_Type[toFill]

## And fill NAs with Unknown
missing <- which(is.na(spp_dat$Revised_fg))
spp_dat$Revised_fg[missing] <- "Unknown"

#################################################
# 9. Sort out units of diversity measures
#################################################

spp_dat$species_Biomassm2 <- NA
levels(spp_dat$WetBiomassUnits)

spp_dat$species_Biomassm2[which(spp_dat$WetBiomassUnits == "g/m2")] <- 
  spp_dat$WetBiomass[which(spp_dat$WetBiomassUnits == "g/m2")]

# mg/m2 -> g/m2 = divide by 1000
spp_dat$species_Biomassm2[which(spp_dat$WetBiomassUnits == "mg/m2")] <- 
  spp_dat$WetBiomass[which(spp_dat$WetBiomassUnits == "mg/m2")] /1000

## convert g to g/m2, divide by the sampled area (in m2)
table(spp_dat$Sampled_Area_Unit[which(spp_dat$WetBiomassUnits == "g")])

spp_dat$species_Biomassm2[which(spp_dat$WetBiomassUnits == "g" & spp_dat$Sampled_Area_Unit == "m2")] <-
  spp_dat$WetBiomass[which(spp_dat$WetBiomassUnits == "g" & spp_dat$Sampled_Area_Unit == "m2")] / spp_dat$Sampled_Area[which(spp_dat$WetBiomassUnits == "g" & spp_dat$Sampled_Area_Unit == "m2")]


## cm2 >- m2 = divide by 10000
spp_dat$species_Biomassm2[which(spp_dat$WetBiomassUnits == "g" & spp_dat$Sampled_Area_Unit == "cm2")] <-
  spp_dat$WetBiomass[which(spp_dat$WetBiomassUnits == "g" & spp_dat$Sampled_Area_Unit == "cm2")] / (spp_dat$Sampled_Area[which(spp_dat$WetBiomassUnits == "g" & spp_dat$Sampled_Area_Unit == "cm2")] / 10000)

summary(spp_dat$species_Biomassm2)
summary(spp_dat$Biomass_fromspecies)
##########################################################
## Abundance values

table(spp_dat$Abundance_Unit)

spp_dat$species_Abundancem2 <- NA
spp_dat$species_Abundancem2[which(spp_dat$Abundance_Unit == "Individuals per m2")] <- spp_dat$Abundance[which(spp_dat$Abundance_Unit == "Individuals per m2")]
## individuals per m3 is basically the same per m2
spp_dat$species_Abundancem2[which(spp_dat$Abundance_Unit == "Individuals per m3")] <- spp_dat$Abundance[which(spp_dat$Abundance_Unit == "Individuals per m3")]


# number of individual when sampled area is in m2
table(spp_dat$Sampled_Area_Unit[which(spp_dat$Abundance_Unit == "Number of individuals")])

spp_dat$species_Abundancem2[which(spp_dat$Abundance_Unit == "Number of individuals" & spp_dat$Sampled_Area_Unit == "m2")] <-
  spp_dat$Abundance[which(spp_dat$Abundance_Unit == "Number of individuals" & spp_dat$Sampled_Area_Unit == "m2")] / spp_dat$Sampled_Area[which(spp_dat$Abundance_Unit == "Number of individuals" & spp_dat$Sampled_Area_Unit == "m2")]


## When in cm2
spp_dat$species_Abundancem2[which(spp_dat$Abundance_Unit == "Number of individuals" & spp_dat$Sampled_Area_Unit == "cm2")] <-
  spp_dat$Abundance[which(spp_dat$Abundance_Unit == "Number of individuals" & spp_dat$Sampled_Area_Unit == "cm2")] / (spp_dat$Sampled_Area[which(spp_dat$Abundance_Unit == "Number of individuals" & spp_dat$Sampled_Area_Unit == "cm2")] / 10000)

# unique(sites$file[which(sites$Site_AbundanceUnits == "Number of individuals" & sites$Sampled_Area_Unit == "m3")])
summary(spp_dat$species_Abundancem2)
#################################################
# 10. Create the dataframe
#################################################

spp_dat$newID <- paste(spp_dat$file.x, spp_dat$Study_Name, spp_dat$Site_Name.x)
sites$newID <- paste(sites$file, sites$Study_Name, sites$Site_Name)

Summary.div <- spp_dat %>% # Start by defining the original dataframe, AND THEN...
  group_by(newID) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    Epi_biomass = sum(species_Biomassm2[Revised_fg == "Epigeic"]),
    Endo_biomass = sum(species_Biomassm2[Revised_fg == "Endogeic"]),
    Ane_biomass = sum(species_Biomassm2[Revised_fg == "Anecic"]),
    EpiEndo_biomass = sum(species_Biomassm2[Revised_fg == "Epi-Endogeic"]),
    Unknown_biomass = sum(species_Biomassm2[Revised_fg == "Unknown"]),
    Epi_abundance = sum(species_Abundancem2[Revised_fg == "Epigeic"]),
    Endo_abundance = sum(species_Abundancem2[Revised_fg == "Endogeic"]),
    Ane_abundance = sum(species_Abundancem2[Revised_fg == "Anecic"]),
    EpiEndo_abundance = sum(species_Abundancem2[Revised_fg == "Epi-Endogeic"]),
    Unknown_abundance = sum(species_Abundancem2[Revised_fg == "Unknown"])
  )

summary.div <- as.data.frame(Summary.div)
str(summary.div)

##########################################################
## Match with site level dataset
##########################################################

sites_fg <- merge(sites, summary.div, by.x = "newID", by.y = "newID", all.x = TRUE)

##########################################################
## Save the data
##########################################################


write.csv(sites_fg, file = file.path(data_out, paste("SiteswithFunctionalGroups_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
