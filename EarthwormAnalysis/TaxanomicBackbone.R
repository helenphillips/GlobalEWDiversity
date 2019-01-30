########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

#################################################
# 1. Loading libraries
#################################################
library(googlesheets)


if(!dir.exists("20_Data")){
  dir.create("20_Data")
}

data_out <- "20_Data"

#################################################
# 2. Load Data
#################################################

# Use the entire dataset that was corrected by Maria and George

data_in_spp <-"9_Data"
files <- list.files(file.path(data_in_spp))
## This is a file I made manually, by combining data from George and Maria (also in this folder)
loadinfg <- "Unique_Species_toSend2018-06-05_MJIB+GB.csv"
spp <- read.csv(file.path(data_in_spp, loadinfg))


## Drilobase can be used for the families, and orders
drilobase <- "DriloBASE_uniqueSpecies"
drilobase <- gs_title(drilobase)
drilo <- as.data.frame(gs_read(drilobase, ws = "Sheet1"))


# Then put into new format

#################################################
# 3. Format Data
#################################################

# Slim down dataframe amd remove duplicates
spp <- spp[,c('Revised', 'Revised_fg', 'Revised_Authority', 'Notes')] ## 304


spp <- spp[!duplicated(spp$Revised), ] # 193
# One row (which is originally three) where the genus and species names are invalid
spp <- spp[!(is.na(spp$Revised)),] # 192

# Functional Groups
levels(spp$Revised_fg)
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "anecic"] <- "Anecic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "endogeic"] <- "Endogeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "endogeic?"] <- "Endogeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "Epi-endogeic"] <- "Epi-Endogeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "epigeic"] <- "Epigeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "N/A"] <- "Unknown"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "endogeic or epiendogeic"] <- "Epi-Endogeic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "Endo-anecic"] <- "Endo-Anecic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "Endogeic, Anecic"] <- "Endo-Anecic"
levels(spp$Revised_fg)[levels(spp$Revised_fg) == "epi-anecic"] <- "Epi-Anecic"

# Genus and species names
spp$Revised <- as.character(spp$Revised)
spp$Genus <-sapply(strsplit(spp$Revised, "\\s+"), "[", 1)
spp$species <-sapply(strsplit(spp$Revised, "\\s+"), "[", 2)

################### 
# Do similar to Drilobase data

head(drilo)
drilo <- drilo[,c('Author of species', 'family', 'genus')]
drilo <- drilo[!duplicated(drilo$genus), ] # 292

#################################################
# 4. Merge data
#################################################

dat <- merge(spp, drilo, by.x= "Genus", by.y = "genus", all.x = TRUE)

dat$phylum <- "Annelida"
dat$class <- "Clitellata"
dat$order <- "Oligochaeta"

#################################################
# 5.  Order data
#################################################

dat <-dat[,c('phylum', 'class', 'order', 'family', 'Genus', 'species', 
             'Revised', 'Revised_fg', 'Author of species', 'Revised_Authority', 'Notes')]


#################################################
# 6. Save data
#################################################

write.csv(dat, file=file.path(data_out, paste("TaxonomicBackbone",  Sys.Date(), ".csv", sep = "")))
