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


########################################################
# 4. Access googledrive
########################################################

x <- gs_ls() ## Authentication


########################################################
# 5. Create dataframe for bibliographic and author info
########################################################

bib_names <- c("Article_Title", "Article_Year", "Article_FirstAuthorSurname",
         "Article_Journal","Article_DOI","DataProvider_Title", "DataProvider_Surname",
         "DataProvider_FirstName", "DataProvider_MiddleInitials", "DataProvider_Email",
         "DataProvider_Institute", "DataProvider_Department", "Number_of_Studies",
         "Total_Number_ofSites","Total_Number_ofSpecies", "Notes","BibKey")

bib <- data.frame(matrix(ncol = length(bib_names), nrow = 0))
colnames(bib) <- bib_names

########################################################
# 6. Get files from googledrive
########################################################

all_files <- x$sheet_title[grep("^\\d*\\_", x$sheet_title, perl = TRUE)]

all_bib <- list()
all_sites <- list()
all_species <- list()

count <- 0

for(file in all_files){
  count <- count + 1
  f <- gs_title(file)
  meta <- as.data.frame(gs_read(f, ws = "MetaData", col_names = FALSE))
  sites <- as.data.frame(gs_read(f, ws = "Site-LevelData"))
  species <- as.data.frame(gs_read(f, ws = "Species-LevelData"))

  ## Sorting out bib dataframe
  meta <- meta[meta[,1] %in% names(bib),] ## To remove unnessecary rows
  meta <- meta[match(names(bib), meta[,1]),] ## To put it in the same order
  bib[file,] <- meta[,2]
  
  all_bib[count] <- bib
  
  ## Adding article ID
  sites <- cbind(file, sites)
  species <- cbind(file, species)
  
  ## Validation of site and species level data
  if(!(all(sites$Study_Name %in% species$Study_ID))) stop("Validation failed: Not all studies have species information")
  if(!(all(sites$Site_Name %in% species$Site_Name))) stop("Validation failed: Not all sites have species information")
  
  sites <- formatSites(sites)
  species <- formatSpecies(species)
  
  ## TODO: Check that dates make sense
  
  
  #### Calculate site level species richness & Check the values there
  juvs <- which(species$LifeStage == "Juvenile")
  notSpecies <- which(is.na(species$SpeciesBinomial) && is.na(species$MorphospeciesID))
  
  if(length(c(juvs, notSpecies)) > 0){
    spR <- as.data.frame(table(species$Study_site[-c(juvs, notSpecies)]))
  } else {spR <- as.data.frame(table(species$Study_site))}
  rm(list=c(juvs, notSpecies))
  names(spR)[2] <- "NumberofSpecies"
  sites <- merge(sites, spR, by.x = "Study_site", by.y = "Var1")
  rm(spR)
  
  ## Calculate site level abundance
  ta <- as.data.frame(tapply(species$Abundance, species$Study_site, sum))
  names(ta) <- "Site_NumberofIndividuals"
  ta$SS <- rownames(ta)
  sites <- merge(sites, ta, by.x = "Study_site", by.y = "SS")
  rm(ta)
  
  ## Calculate site level biomass
  ## TODO
  
  ## Check if site level species richness values were given
  check <- which(!(is.na(sites$SpeciesRichness)) && (sites$SpeciesRichnessUnit == "Number of species"))
  if(length(check) > 0){
    if(any(sites$SpeciesRichness[check] != sites$NumberofSpecies[check])){cat(paste("\n", file, ":Some of the site level species richness values do not add up"))}
  }
  

  ## Now to make a species level dataframe with all the variables in
  site_species <- merge(species, sites, by.x = "Study_site", by.y = "Study_site", all.x = TRUE)
  
  all_sites[count] <- sites
  all_species[count] <- site_species
  
}
