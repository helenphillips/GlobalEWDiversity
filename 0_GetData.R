########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

########################################################
# 2. Create folder if it doesn't exist to save data into
########################################################

if(!dir.exists("0_Data")){
  dir.create("0_Data")
}
data_out <- "0_Data"

########################################################
# 3. Libraries
########################################################

library(googlesheets)

source("Functions/FormatData.R")
source("Functions/CalculateSitelevelMetrics.R")

########################################################
# 4. Access googledrive
########################################################

x <- gs_ls() ## Authentication



########################################################
# 5. Create dataframe for bibliographic and author info
########################################################
template <- "DatasetTemplateoriginal"
template_meta <- gs_title(template)
meta <- as.data.frame(gs_read(template_meta, ws = "MetaData", col_names = FALSE))

bib_names <- meta[,1]

bib <- data.frame(matrix(ncol = length(bib_names) + 1, nrow = 0)) ## plus one so I can add a file name
colnames(bib) <- c(bib_names, "file")



######################################################
# 5.5. Create template of site info for checking
#######################################################

sitetemplate <- as.data.frame(gs_read(template_meta, ws = "Site-LevelData"))
speciestemplate <- as.data.frame(gs_read(template_meta, ws = "Species-LevelData"))


########################################################
# 6. Get files from googledrive
########################################################

all_files <- x$sheet_title[grep("^\\d*\\_", x$sheet_title, perl = TRUE)]

cat(paste("\nFound", length(all_files), "datasheets"))

options( warn = 2 )

#all_bib <- list(length = length(bib_names))
#names(all_bib) <- bib_names
all_sites <- list()
all_species <- list()

########################################################
# 7. Start processing data
########################################################

count <- 0

for(file in all_files){
  count <- count + 1
  f <- gs_title(file)
  meta <- as.data.frame(gs_read(f, ws = "MetaData", col_names = FALSE))
  sites <- as.data.frame(gs_read(f, ws = "Site-LevelData"))
  species <- as.data.frame(gs_read(f, ws = "Species-LevelData"))
  
  ## Add file info to metadata
  meta[nrow(meta) + 1,] <- c("file", file)
  
  
  ### Check that most recent version of template, else skip it
  if(!(all(names(bib) %in% meta[,1]))) {
    print(paste(file, ": failed on meta sheet"))
    next }
  if(!(all(names(sitetemplate) %in% names(sites)))){
    print(paste(file, ": failed on sites sheet"))
    next }
  if(!(all(names(speciestemplate) %in% names(species)))) {
    print(paste(file, ": failed on species sheet"))
    next }
  
  ## Sorting out bib dataframe
  meta <- meta[meta[,1] %in% names(bib),] ## To remove unnessecary rows
  meta <- meta[match(names(bib), meta[,1]),] ## To put it in the same order
  bib[count,] <- meta[,2]

  
  ## Adding article ID
  sites <- cbind(file, sites)
  if(nrow(species) > 0){species <- cbind(file, species)}
  
  ## Formatting
  sites <- formatSites(sites)
  if(nrow(species) > 0){species <- formatSpecies(species)}
  
  
  ## TODO: Check that dates make sense
  
  
  #### Calculate site level species richness from species list & Check the values there
  if(nrow(species) > 0){
    juvs <- which(species$LifeStage == "Juvenile")
    notSpecies <- which(is.na(species$SpeciesBinomial) & is.na(species$MorphospeciesID))
  
    if(length(c(juvs, notSpecies)) > 0){
      spR <- as.data.frame(table(species$Study_site[-c(juvs, notSpecies)]))
    } else {spR <- as.data.frame(table(species$Study_site))}
    rm(juvs)
    rm(notSpecies)
    names(spR)[2] <- "NumberofSpecies"
    sites <- merge(sites, spR, by.x = "Study_site", by.y = "Var1")
    rm(spR)
  
  ## Calculate site level abundance
    ta <- as.data.frame(tapply(species$Abundance, species$Study_site, sum))
    names(ta) <- "Individuals_fromspecies"
    ta$Individuals_fromspeciesUnits <- species$Abundance_Unit[1]
    ta$SS <- rownames(ta)
    sites <- merge(sites, ta, by.x = "Study_site", by.y = "SS")
    rm(ta)
  
  ## Calculate site level biomass
    bm <- as.data.frame(tapply(species$WetBiomass, species$Study_site, sum))
    names(bm) <- "Biomass_fromspecies"
    bm$Biomass_fromspeciesUnits <- species$WetBiomassUnits[1]
    bm$SS <- rownames(bm)
    sites <- merge(sites, bm, by.x = "Study_site", by.y = "SS")
    rm(bm)
    
  ## Check if site level species richness values were given
    check <- which(!(is.na(sites$SpeciesRichness)) && (sites$SpeciesRichnessUnit == "Number of species"))
    if(length(check) > 0){
      if(any(sites$SpeciesRichness[check] != sites$NumberofSpecies[check])){cat(paste("\n", file, ":Some of the site level species richness values do not add up\n"))}
    }
    rm(check)
    
  ## Check if abundance and biomass values were given
    check <- which(!(is.na(sites$Site_WetBiomass)) & !(is.na(sites$Biomass_fromspecies)))
    if(length(check) > 0){
      if(any(sites$Site_WetBiomass[check] != sites$Biomass_fromspecies[check])){cat(paste("\n", file, ":Some of the site level biomass values do not add up\n"))}
    }
    rm(check)
    
    check <- which(!(is.na(sites$Site_Abundance)) & !(is.na(sites$Individuals_fromspecies)))
    if(length(check) > 0){
      if(any(sites$Site_Abundance[check] != sites$Individuals_fromspecies[check])){cat(paste("\n", file, ":Some of the site level abundance values do not add up\n"))}
    }
    rm(check)
    
    
  }else{
    sites$NumberofSpecies <- NA
    sites$Individuals_fromspecies <- NA
    sites$Individuals_fromspeciesUnits <- NA
    sites$Biomass_fromspecies <- NA
    sites$Biomass_fromspeciesUnits <- NA}
  
  ## Now to make a species level dataframe with all the variables in
  if(nrow(species) > 0){site_species <- merge(species, sites, by.x = "Study_site", by.y = "Study_site", all.x = TRUE)
  all_species[[count]] <- site_species
    }
  
  
  all_sites[[count]] <- sites
  
  
}

########################################################
# 8. Create two complete dataframes
########################################################


sites <- do.call("rbind", all_sites)
species <- do.call("rbind", all_species)


########################################################
# 9. Amalgamte site level values of species richness, biomass and abundance
########################################################

sitelevel_spR <- which(is.na(sites$SpeciesRichness) & !(is.na(sites$NumberofSpecies)))
sites$SpeciesRichness[sitelevel_spR] <-sites$NumberofSpecies[sitelevel_spR]
sites$SpeciesRichnessUnit[sitelevel_spR] <- "Number of species"
rm(sitelevel_spR)

sites$Site_WetBiomassUnits <- as.character(sites$Site_WetBiomassUnits)

sitelevel_biom <- which(is.na(sites$Site_WetBiomass) & !(is.na(sites$Biomass_fromspecies)))
if(length(sitelevel_biom ) > 0) {
  sites$Site_WetBiomass[sitelevel_biom] <-sites$Biomass_fromspecies[sitelevel_biom]
  sites$Site_WetBiomassUnits[sitelevel_biom] <- as.character(sites$Biomass_fromspeciesUnits[sitelevel_biom])
  }
rm(sitelevel_biom)

sites$Site_AbundanceUnits <- as.character(sites$Site_AbundanceUnits)


sitelevel_abund <- which(is.na(sites$Site_Abundance) & !(is.na(sites$Individuals_fromspecies)))
if(length(sitelevel_abund ) > 0) {
  sites$Site_Abundance[sitelevel_abund] <-sites$Individuals_fromspecies[sitelevel_abund]
  sites$Site_AbundanceUnits[sitelevel_abund] <- as.character(sites$Individuals_fromspeciesUnits[sitelevel_abund])
  }
rm(sitelevel_abund)


colsToRemove <- c("NumberofSpecies", "Individuals_fromspecies", "Individuals_fromspeciesUnits", "Biomass_fromspecies", "Biomass_fromspeciesUnits")
sites[,names(sites) %in% colsToRemove] <- NULL


########################################################
# 10. Save the data
########################################################
sites$ID <- 1:nrow(sites)
sites$ID <- paste(sites$file, sites$Study_Name, sites$ID, sep = "_")

write.csv(sites, file = file.path(data_out, paste("sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
write.csv(species, file = file.path(data_out, paste("species_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
write.csv(bib, file = file.path(data_out, paste("Metadata_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

########################################################
# 11. Save the data for Carlos
########################################################

# carlos_sites <- sites[,c("ID", "Latitude__decimal_degrees", "Longitude__Decimal_Degrees")]
# names(carlos_sites) <- c("ID", "LAT", "LONG")
# write.csv(carlos_sites, file = file.path(data_out, paste("sitesIDLatLong_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)


########################################################
# 11. Save the data for Erin and Olga
########################################################

# eo_sites <- sites[1:200,]
# 
# eo_species <- species[7500:8159,]
# 
# toDel <- c("file.y", "Study_Name", "Site_Name.y", "SpeciesRichness", 
#            "SpeciesRichnessUnit","Site_WetBiomass","Site_WetBiomassUnits",   
#            "Site_Abundance","Site_AbundanceUnits","NumberofSpecies",
#            "Individuals_fromspecies","Individuals_fromspeciesUnits", 
#            "Biomass_fromspecies","Biomass_fromspeciesUnits")
# 
# eo_species[,which(names(eo_species) %in% toDel)] <- NULL
# 
# write.csv(eo_species, file = file.path(data_out, paste("speciesExample_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
# write.csv(eo_sites, file = file.path(data_out, paste("sitesExample", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
