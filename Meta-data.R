if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}
#################################################
# 1. Libraries
#################################################
library(lme4)

#################################################
# 3. Create directories
#################################################


#################################################
# 4. Load in models
#################################################

models <- "Models"


load(file.path(models, "richnessmodel.rds"))
load(file.path(models, "biomassmodel_full.rds"))
load(file.path(models, "abundancemodel_full.rds"))
# load(file.path(models, "fgrichnessmodel.rds"))

data_in <- "4_Data"


#################################################
# 4. Load in data
#################################################
files <- list.files(file.path(data_in))
files <- files[grep("sitesRichness_2", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

richness <- read.csv(file.path(data_in, loadin))

files <- list.files(file.path(data_in))
files <- files[grep("sitesBiomass_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

biomass <- read.csv(file.path(data_in, loadin))

files <- list.files(file.path(data_in))
files <- files[grep("sitesAbundance_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

abundance <- read.csv(file.path(data_in, loadin))

#################################################
# 4. Load in data
#################################################
names(richness)[names(richness) == "scaleElevation"] <- "ScaleElevation"

alldat <- rbind(richness, biomass, abundance)

## unique sites
length(unique(alldat$ID))
## 7051
length(unique(alldat$Study_site))
## 7048


length(unique(alldat$Study_Name)) # 229
length(unique(alldat$file)) # 181


countries <- unique(alldat$Country)
## The NAs are USA and spain...both already represented

countries[countries == "MEXICO"] <- "Mexico"
countries[countries == "United States"] <- "USA"
countries[countries == "United Kingdom"] <- "UK"
countries[countries == "Pueto Rico"] <- "Puerto Rico"
countries[countries == "Hawaii"] <- "USA"
countries[countries == "Wales"] <- "UK"

countries <- unique(countries)

countries <- droplevels(as.factor(countries))
levels(countries)





###############################################
## Load in the bib data
################################################

bibdata <- "0_Data"

files <- list.files(file.path(bibdata))
files <- files[grep("Metadata_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

bib <- read.csv(file.path(bibdata, loadin))


files <- unique(alldat$file)

bib <- droplevels(bib[bib$file %in% files,])

table(bib$Data.From.Paper)


######## published or not?
one <- bib[grep("unpublished", bib$Article_Title, ignore.case = TRUE),]
two <- bib[grep("thesis", bib$Article_Title, ignore.case = TRUE),]
three <- bib[is.na(bib$Article_Title),]
unique(bib$Article_Year) ## there's no unpublished in the article year

t <- rbind(one, two, three)

length(unique(t$file))


######################################################
## A UNIQUE DATASET
########################################################
uniqueDat <- alldat
alldat[,grep("scale", names(uniqueDat), ignore.case = TRUE)] <- NULL
uniqueDat <- uniqueDat[!duplicated(uniqueDat$Study_site), ]
## The scaled variables would be different for each dataframe, so removing them
nrow(uniqueDat)


summary(uniqueDat$PH)
nrow(uniqueDat) - 2628
summary(uniqueDat$CEC)
nrow(uniqueDat) -  6546


summary(uniqueDat$Base_Saturation_percent)
summary(uniqueDat$Organic_Carbon__percent)
summary(uniqueDat$Soil_Organic_Matter__percent)
summary(uniqueDat$C.N_ratio)
summary(uniqueDat$Sand__percent)
summary(uniqueDat$Silt__percent)
summary(uniqueDat$Clay__percent)
summary(uniqueDat$USDA_SoilTexture)
summary(uniqueDat$Soil_Moisture_percent)
summary(uniqueDat$WRB_FAO_SoilType)

summary(uniqueDat$Soil_Organic_Matter__percent)

#####################
# Studies with soil
####################

library(plyr)
library(dplyr)

########
##Studies with identical climate variables
#######

soil.df <- uniqueDat %>% # Start by defining the original dataframe, AND THEN...
  group_by(Study_Name) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    PH = mean(PH, na.rm = TRUE),
    CEC = mean(CEC, na.rm = TRUE),
    Organic_Carbon__percent = mean(Organic_Carbon__percent, na.rm = TRUE),
    Soil_Organic_Matter__percent = mean(Soil_Organic_Matter__percent, na.rm = TRUE),
    C.N_ratio = mean(C.N_ratio, na.rm = TRUE),
    Sand__percent = mean(Sand__percent, na.rm = TRUE)
  )

soildf <- as.data.frame(soil.df)

someSoil <- soildf[apply(soildf[c('PH','CEC','Organic_Carbon__percent', 'Soil_Organic_Matter__percent', 'C.N_ratio', 'Sand__percent')],1,function(x) any(x > 0)),]
length(unique(someSoil$Study_Name)) ## 182 with one NA
#######################
## Studies with no variation in climate
########################

var.df <- uniqueDat %>% # Start by defining the original dataframe, AND THEN...
  group_by(Study_Name) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    bio10_1 = var(bio10_1),
    bio10_4 = var(bio10_4),
    bio10_7 = var(bio10_7),
    bio10_12 = var(bio10_12),
    bio10_15 = var(bio10_15)
  )

vardf <- as.data.frame(var.df)

someClimate <- vardf[apply(vardf [c('bio10_1','bio10_4','bio10_7', 'bio10_12', 'bio10_15')],1,function(x) any(x > 0)),]
noClimate <- vardf[apply(vardf [c('bio10_1','bio10_4','bio10_7', 'bio10_12', 'bio10_15')],1,function(x) all(x == 0)),]


#########################################
## studies with high diversity
#########################################


highRichness <- droplevels(richness[richness$SpeciesRichness > 10,])
unique(highRichness$Study_Name)


highAbundance <- droplevels(abundance[abundance$Sites_Abundancem2 > 600,])
unique(highAbundance$Study_Name)


highBiomass <- droplevels(biomass[biomass$Site_WetBiomass > 300,])
unique(highBiomass$Study_Name)
