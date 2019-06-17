
########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

source("Functions/FormatData.R")
########################################################
# 2. Create folder if it doesn't exist to save data into
########################################################

if(!dir.exists("5_Data")){
  dir.create("5_Data")
}
data_out <- "5_Data"

#################################################
# 3. Loading in variables
#################################################

data_in <-"3_Data"
files <- list.files(file.path(data_in))
files <- files[grep("sitesWithChelsaAndSoilAndOthers_", files)]


file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinsites <- loadin[grep("sitesWithChelsaAndSoilAndOthers_", loadin)]

bib_in <-"0_Data"
files <- list.files(file.path(bib_in))
files <- files[grep("Metadata_", files)]
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
# 6. Biomass and Abundance units
#################################################
sites$Site_Biomassm2 <- NA
sites$Site_Biomassm2[which(sites$Site_WetBiomassUnits == "g/m2")] <- sites$Site_WetBiomass[which(sites$Site_WetBiomassUnits == "g/m2")]


# mg/m2 -> g/m2 = divide by 1000
table(sites$Site_WetBiomassUnits)
sites$Site_Biomassm2[which(sites$Site_WetBiomassUnits == "mg/m2")] <- sites$Site_WetBiomass[which(sites$Site_WetBiomassUnits == "mg/m2")] /1000


## convert g to g/m2, divide by the sampled area (in m2)


table(sites$Sampled_Area_Unit[which(sites$Site_WetBiomassUnits == "g")])

sites$Site_Biomassm2[which(sites$Site_WetBiomassUnits == "g" & sites$Sampled_Area_Unit == "m2")] <-
  sites$Site_WetBiomass[which(sites$Site_WetBiomassUnits == "g" & sites$Sampled_Area_Unit == "m2")] / sites$Sampled_Area[which(sites$Site_WetBiomassUnits == "g" & sites$Sampled_Area_Unit == "m2")]


## cm2 >- m2 = divide by 10000
table(sites$Sampled_Area_Unit[which(sites$Site_WetBiomassUnits == "g")])

sites$Site_Biomassm2[which(sites$Site_WetBiomassUnits == "g" & sites$Sampled_Area_Unit == "cm2")] <-
  sites$Site_WetBiomass[which(sites$Site_WetBiomassUnits == "g" & sites$Sampled_Area_Unit == "cm2")] / (sites$Sampled_Area[which(sites$Site_WetBiomassUnits == "g" & sites$Sampled_Area_Unit == "cm2")] / 10000)

summary(sites$Site_Biomassm2)
summary(sites$Site_WetBiomass)
##########################################################
## Abundance values

table(sites$Site_AbundanceUnits)

sites$Sites_Abundancem2 <- NA
sites$Sites_Abundancem2[which(sites$Site_AbundanceUnits == "Individuals per m2")] <- sites$Site_Abundance[which(sites$Site_AbundanceUnits == "Individuals per m2")]
## individuals per m3 is basically the same per m2
sites$Sites_Abundancem2[which(sites$Site_AbundanceUnits == "Individuals per m3")] <- sites$Site_Abundance[which(sites$Site_AbundanceUnits == "Individuals per m3")]


# number of individual when sampled area is in m2
table(sites$Sampled_Area_Unit[which(sites$Site_AbundanceUnits == "Number of individuals")])

sites$Sites_Abundancem2[which(sites$Site_AbundanceUnits == "Number of individuals" & sites$Sampled_Area_Unit == "m2")] <-
  sites$Site_Abundance[which(sites$Site_AbundanceUnits == "Number of individuals" & sites$Sampled_Area_Unit == "m2")] / sites$Sampled_Area[which(sites$Site_AbundanceUnits == "Number of individuals" & sites$Sampled_Area_Unit == "m2")]


## When in cm2
sites$Sites_Abundancem2[which(sites$Site_AbundanceUnits == "Number of individuals" & sites$Sampled_Area_Unit == "cm2")] <-
  sites$Site_Abundance[which(sites$Site_AbundanceUnits == "Number of individuals" & sites$Sampled_Area_Unit == "cm2")] / (sites$Sampled_Area[which(sites$Site_AbundanceUnits == "Number of individuals" & sites$Sampled_Area_Unit == "cm2")] / 10000)

# unique(sites$file[which(sites$Site_AbundanceUnits == "Number of individuals" & sites$Sampled_Area_Unit == "m3")])

hist(sites$Site_Biomassm2)
hist(sites$Sites_Abundancem2)
summary(sites$Sites_Abundancem2)
summary(sites$Site_Abundance)

# tail(sites[which(is.na(sites$Sites_Abundancem2) & !(is.na(sites$Site_Abundance))),c(2, 15, 16, 59, 60, 96)], 90)

##### One study has ridiculous number of individuals per m2.
## I have checked it, and I have transferred correctly from paper.
## Removing, as massive outlier.

sites$file[which(sites$Sites_Abundancem2 > 3000)]

sites <- droplevels(sites[sites$file != "949_Johnson-Maynard2002",]) # 7805
hist(sites$Site_Biomassm2)
hist(sites$Sites_Abundancem2)
summary(sites$Site_Biomassm2)
summary(sites$Sites_Abundancem2)

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

# nrow(sites[!(is.na(sites$Site_WetBiomass)) && is.na(sites$Site_WetBiomassUnits),])

any(!(is.na(sites$Site_Abundance)) && is.na(sites$Site_AbundanceUnits))



#################################################
#  Save file
#################################################

write.csv(sites, file = file.path(data_out, paste("sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

## Sites for Carlos
# write.csv(sites, file = file.path(data_out, paste("sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)

