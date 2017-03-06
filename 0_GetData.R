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
  
  
  ## TODO: Check that dates make sense
  
  
  ## Add in "Biomass" column if not present
  if(!("WetBiomass" %in% names(sites))){
    sites['WetBiomass'] <- NA
    sites['WetBiomassUnits'] <- NA}
  if(!("WetBiomass" %in% names(species))){
    species['WetBiomass'] <- NA
    species['WetBiomassUnits'] <- NA}
 
  
  #### Calculate site level species richness & Check the values there
  table(species$Site_Name[which(species$LifeStage != "Juvenile")])
  
  #### Calculate site level abundance & Check the values there
  
  
  #### Calculate site level biomass & Check the values there
  
  
  
  all_sites[count] <- sites
  all_species[count] <- species
  
}
