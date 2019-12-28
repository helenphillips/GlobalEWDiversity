
########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}


if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\Users\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
  
}


source("Functions/FormatData.R")

library(dplyr)
library(reshape2)
library(plyr)
########################################################
# 2. Create folder if it doesn't exist to save data into
########################################################

if(!dir.exists("16_Data")){
  dir.create("16_Data")
}
data_out <- "16_Data"


#################################################
# 3. Loading in variables
#################################################
## Species data

data_in_spp <-"15_Data"
files <- list.files(file.path(data_in_spp))
## This is a file I made manually, by combining data from George and Maria (also in this folder)
loadinfg <- "UniqueSpeciestoSend_2019-06-19_Final.csv"
## Some species identifications are in press - look at George's final file to identify. 
# Guillaume Rosseau would need to be contacted to get the ecological categories


## Site data


data_in_sites <-"7_Data"

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
spp_dat <- read.csv(file.path(data_in_spp2, loadin_spp)) # 21210

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
## revised species names to append to the original

#################################################
# 7. Merge with species level dataset
#################################################

spp_dat <- merge(spp_dat, spp, by.x = "SpeciesBinomial", by.y = "SpeciesBinomial", all.x = TRUE) # 21210
## Merge using the old names

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
# 9.5 Save this dataset
#################################################

write.csv(spp_dat, file = file.path(data_out, paste("SpecieswithFunctionalGroups_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

#################################################
# 10. Create the dataframe
#################################################

spp_dat$newID <- paste(spp_dat$file.x, spp_dat$Study_Name, spp_dat$Site_Name.x)
sites$newID <- paste(sites$file, sites$Study_Name, sites$Site_Name)
detach(package:plyr) 
Summary.div <- spp_dat %>% # Start by defining the original dataframe, AND THEN...
  group_by(newID) %>% # Define the grouping variable, AND THEN...
  summarize( # Now you define your summary variables with a name and a function...
    Epi_biomass = sum(species_Biomassm2[which(Revised_fg == "Epigeic")], na.rm = TRUE),
    Endo_biomass = sum(species_Biomassm2[which(Revised_fg == "Endogeic")], na.rm = TRUE),
    Ane_biomass = sum(species_Biomassm2[which(Revised_fg == "Anecic")], na.rm = TRUE),
    EpiEndo_biomass = sum(species_Biomassm2[which(Revised_fg == "Epi-Endogeic")], na.rm = TRUE),
    Unknown_biomass = sum(species_Biomassm2[which(Revised_fg == "Unknown")], na.rm = TRUE),
    Epi_abundance = sum(species_Abundancem2[which(Revised_fg == "Epigeic")], na.rm = TRUE),
    Endo_abundance = sum(species_Abundancem2[which(Revised_fg == "Endogeic")], na.rm = TRUE),
    Ane_abundance = sum(species_Abundancem2[which(Revised_fg == "Anecic")], na.rm = TRUE),
    EpiEndo_abundance = sum(species_Abundancem2[which(Revised_fg == "Epi-Endogeic")], na.rm = TRUE),
    Unknown_abundance = sum(species_Abundancem2[which(Revised_fg == "Unknown")], na.rm = TRUE)
  )

summary.div <- as.data.frame(Summary.div)
str(summary.div)


#######################################################
## SPECIES RICHNESS
########################################################

juvs <- which(spp_dat$LifeStage == "Juvenile")
notSpecies <- which(is.na(spp_dat$Revised) & is.na(spp_dat$MorphospeciesID))

notSp <- union(juvs, notSpecies)
spR <- spp_dat[-notSp,]
library(plyr)
t <- ddply(spR, c("newID", "Revised_fg"), summarise, nrows = length(Revised_fg))
t2 <- dcast(t, newID ~ Revised_fg, value.var = "nrows")
names(t2)[2:6] <- c("Ane_richness", "Endo_richness", "EpiEndo_richness", "Epi_richness", "Unknown_richness")
t2$total <- rowSums(t2[,2:ncol(t2)], na.rm = TRUE)
min(t2$total) 

## so The NAs should be 0s

t2[is.na(t2)] <- 0
t2$total <- NULL
###################
## Calculate the functional group richness as well
###################

t3 <- t2
t3$Unknown_richness <- NULL
t3$Ane_richness <- ifelse(t2$Ane_richness == 0, 0, 1)
t3$Endo_richness <- ifelse(t2$Endo_richness == 0, 0, 1)
t3$EpiEndo_richness <- ifelse(t2$EpiEndo_richness == 0, 0, 1)
t3$Epi_richness <- ifelse(t2$Epi_richness == 0, 0, 1)

t3$FGRichness <- rowSums(t3[,2:5])
t3 <- t3[,c('newID', 'FGRichness')]
t2 <- merge(t2, t3, by = "newID")

##########################################################
## Match with site level dataset
##########################################################

sites_fg <- merge(sites, summary.div, by.x = "newID", by.y = "newID", all.x = TRUE)
sites_fg <- merge(sites_fg, t2, by.x = "newID", by.y = "newID", all.x = TRUE)


##########################################################
## Fill in missing data
##########################################################
# There may be studies where some functional groups were present
# but not others (or sites nothing at all was found, within a study)
# Need to ensure that they are marked as zeros, not NAs

# I am tired, so let's do it a safe way


sites_fg2 <- c()

all_studies <- unique(sites_fg$studyID)

for(s in 1:length(all_studies)){
  
  print(all_studies[s])
  
  temp <- sites_fg[sites_fg$studyID == all_studies[s],]
  
  vars <- c("biomass", "abundance", "_richness")
  
  for(v in 1:length(vars)){
    print(vars[v])
    bio <- grep(vars[v], names(temp))
    
    if(any(!(is.na(temp[,bio])))){ # If any of them are not NA
      # Then make no sure no NAs
      print("There are some values")
      if(any(is.na(temp[,bio]))){ # If there are NAs
        print("There are some NAs")
        for( col in 1:length(bio)){
          print("removing NAs")
          temp[,bio[col]] <- ifelse(is.na(temp[,bio[col]]), 0, temp[,bio[col]])
        }
        
        
      }
      
    }
  }
  
  sites_fg2 <- rbind(sites_fg2, temp)
  
}




##########################################################
## Save the data
##########################################################


write.csv(sites_fg2, file = file.path(data_out, paste("SiteswithFunctionalGroups_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

