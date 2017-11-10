
########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

source("Functions/FormatData.R")
########################################################
# 2. Create folder if it doesn't exist to save data into
########################################################

if(!dir.exists("2_Data")){
  dir.create("2_Data")
}
data_out <- "2_Data"

#################################################
# 3. Loading in variables
#################################################

data_in <-"1_Data"
files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinsites <- loadin[grep("sites_", loadin)]

bib_in <-"0_Data"
files <- list.files(file.path(bib_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinbib <- loadin[grep("Metadata_", loadin)]


#################################################
# 4. Load in data
#################################################

sites <- read.csv(file.path(data_in, loadinsites))
bib <- read.csv(file.path(bib_in, loadinbib))

#################################################
# 5. Get rid of studies with selected species 
#################################################
bib$Entire.Community <- as.factor(bib$Entire.Community)
all_spp <- bib$file[which(bib$Entire.Community != "no - select species sampled")]
all_spp <- c(as.vector(all_spp), as.vector(bib$file[which(is.na(bib$Entire.Community))]))
sites <- sites[sites$file %in% all_spp,]

#################################################
# 6. MAke some factor levels consistent
#################################################


#################################################
# 7. Rename factor levels
#################################################
sites$Management_System <- as.character(sites$Management_System)
sites$Management_System[which(is.na(sites$Management_System))] <- "None"
sites$Management_System <- as.factor(sites$Management_System)

sites$LandUse <- as.character(sites$LandUse)
sites$LandUse[which(is.na(sites$LandUse))] <- "Unknown"
sites$LandUse <- as.factor(sites$LandUse)

sites$HabitatCover <- as.character(sites$HabitatCover)
sites$HabitatCover[which(is.na(sites$HabitatCover))] <- "Unknown"
sites$HabitatCover <- as.factor(sites$HabitatCover)



#################################################
# 8.0.1 Create a new variable
#################################################
keep <- c("Primary vegetation","Secondary vegetation","Urban","Unknown")
sites$LU_Mgmt <- as.factor(ifelse(sites$LandUse %in% keep, as.character(sites$LandUse), as.character(sites$Management_System)))


#################################################
# 8.0.2 Create a new ESA variable
#################################################
sites$ESA <- sites$HabitatCover

sites$ESA <- as.character(sites$ESA)



prod_herb <- which(sites$LandUse == "Production - Arable" | sites$Management_System == "Annual crop")
sites$ESA[prod_herb] <- "Production - Herbaceous"

prod_planta <- which(sites$LandUse == "Production - Crop plantations" | sites$LandUse == "Production - Wood plantation" | sites$Management_System == "Tree plantations")
sites$ESA[prod_planta] <- "Production - Plantation"

integratessys <- which(sites$Management_System == "Integrated systems")
sites$ESA[integratessys] <- "Cropland/Other vegetation mosaic"

pastures <- which(sites$Management_System == "Pastures (grazed lands)")
sites$ESA[pastures] <- "Herbaceous"

# unique(sites$ESA[sites$LandUse == "Pasture"])

table(sites$ESA)



#### There are some empty cells
## And at the moment, can't do anything with them
# nodata <- droplevels(sites[is.na(sites$ESA),])

sites$ESA[is.na(sites$ESA)] <- "Unknown"





#################################################
# 8. Set reference levels
#################################################

sites <- SiteLevels(sites) 


#################################################
# 8.1 Check that all land uses have comparisons between studies
#################################################
# Which land uses do we have comparisons of within a study
known<- sites[sites$LandUse != "Unknown",]
morethan1 <- names(which(tapply(known$LandUse, known$file, function(x) length(unique(x))) > 1))

landusecomp <- unique(known$LandUse[known$file %in% morethan1]) ## Any land uses in this list, will have a comparison
fileswith <- unique(known$file[!(any(known$LandUse %in% landusecomp))]) ## which sources have sites which do not have a comparison


#################################################
# 8.2 Check if any habitat covers have comaprisons
#################################################
# Which land uses do we have comparisons of within a study
known<- sites[sites$LandUse != "Unknown",]
morethan1 <- names(which(tapply(known$LandUse, known$file, function(x) length(unique(x))) > 1))

landusecomp <- unique(known$LandUse[known$file %in% morethan1]) ## Any land uses in this list, will have a comparison
fileswith <- unique(known$file[!(any(known$LandUse %in% landusecomp))]) ## which sources have sites which do not have a comparison

#################################################
# 8.3 Check if all Land-use/Management categories have comaprisons
#################################################
known<- sites[sites$LU_Mgmt != "Unknown",]
morethan1 <- names(which(tapply(known$LU_Mgmt, known$file, function(x) length(unique(x))) > 1))

landusecomp <- unique(known$LU_Mgmt[known$file %in% morethan1]) ## Any land uses in this list, will have a comparison
fileswith <- unique(known$file[!(any(known$LU_Mgmt %in% landusecomp))]) ## which sources have sites which do not have a comparison

rm(known);rm(morethan1); rm(landusecomp); rm(fileswith)
#################################################
# 9. Check that all biomass and abundance values have units
#################################################

any(!(is.na(sites$Site_WetBiomass)) && is.na(sites$Site_WetBiomassUnits)) ## If true, there's no units for samples
any(!(is.na(sites$Site_Abundance)) && is.na(sites$Site_AbundanceUnits))


#################################################
#  Save file
#################################################

write.csv(sites, file = file.path(data_out, paste("sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

## Sites for Carlos
write.csv(sites, file = file.path(data_out, paste("sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

