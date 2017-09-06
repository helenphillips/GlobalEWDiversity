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

bib <- data.frame(matrix(ncol = length(bib_names), nrow = 0))
colnames(bib) <- bib_names



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

#all_bib <- list(length = length(bib_names))
#names(all_bib) <- bib_names
all_sites <- list()
all_species <- list()


count <- 0

for(file in all_files){
  count <- count + 1
  f <- gs_title(file)
  meta <- as.data.frame(gs_read(f, ws = "MetaData", col_names = FALSE))
  sites <- as.data.frame(gs_read(f, ws = "Site-LevelData"))
  species <- as.data.frame(gs_read(f, ws = "Species-LevelData"))
  
  ### Check that most recent version of template, else skip it
  if(!(all(names(bib) %in% meta[,1]))) next
  if(!(all(names(sitetemplate) %in% names(sites)))) next
  if(!(all(names(speciestemplate) %in% names(species)))) next
  
  ## Sorting out bib dataframe
  meta <- meta[meta[,1] %in% names(bib),] ## To remove unnessecary rows
  meta <- meta[match(names(bib), meta[,1]),] ## To put it in the same order
  bib[1,] <- meta[,2]
  #all_bib[[count]] <- bib[1,]
  
  ## Adding article ID
  sites <- cbind(file, sites)
  if(nrow(species) > 0){species <- cbind(file, species)}
  
  
  sites <- formatSites(sites)
  if(nrow(species) > 0){species <- formatSpecies(species)}
  
  ## Validation of site and species level data
  # already_checked <- c("4643_Nieminen2011", "230_Falco2015", "279_Pansu2015", "665_Raty2004", "000_Guernion2014")
  # if(!(file %in% already_checked)){
  #  if(!(all(levels(sites$Study_site) %in% levels(species$Study_site)))) stop("Validation failed: Not all sites have species information")
  # }
  # rm(already_checked)
  
  
  ## TODO: Check that dates make sense
  
  
  #### Calculate site level species richness from species list & Check the values there
  if(nrow(species) > 0){
    juvs <- which(species$LifeStage == "Juvenile")
    notSpecies <- which(is.na(species$SpeciesBinomial) && is.na(species$MorphospeciesID))
  
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
      if(any(sites$SpeciesRichness[check] != sites$NumberofSpecies[check])){cat(paste("\n", file, ":Some of the site level species richness values do not add up"))}
    }
    rm(check)
    
  ## Check if abundance and biomass values were given
    
    
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



sites <- do.call("rbind", all_sites)
species <- do.call("rbind", all_species)

write.csv(sites, file = file.path(data_out, paste("sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
write.csv(species, file = file.path(data_out, paste("species_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
