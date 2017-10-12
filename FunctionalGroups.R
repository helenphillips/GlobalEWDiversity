
library(dplyr)


##########################################################
## Load in data
##########################################################

data_in <-"0_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
loadin <- files[grep(date, files)]
loadinspecies <- loadin[grep("species_", loadin)]

species <- read.csv(file.path(data_in, loadinspecies))



data_in <-"3_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
loadin <- files[grep(date, files)]

sites <- read.csv(file.path(data_in, loadin))

##########################################################
## Calculate site level measures of the three main functional groups
##########################################################


Summary.df <- species %>% # Start by defining the original dataframe, AND THEN...
  group_by(Study_site) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    Epi_biomass = sum(WetBiomass[Functional_Type == "Epigeic"]),
    Epi_abundance = sum(Abundance[Functional_Type == "Epigeic"]),
    Endo_biomass = sum(WetBiomass[Functional_Type == "Endogeic"]),
    Endo_abundance = sum(Abundance[Functional_Type == "Endogeic"]),
    Ane_biomass = sum(WetBiomass[Functional_Type == "Anecic"]),
    Ane_abundance = sum(Abundance[Functional_Type == "Anecic"])
  )

summary.df <- as.data.frame(Summary.df)
summary.df <- summary.df[complete.cases(summary.df),]

##########################################################
## Match with site level dataset
##########################################################

sites_fg <- (merge(sites, summary.df, by = "Study_site", all.x = TRUE))


##########################################################
## simple model
##########################################################

