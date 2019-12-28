if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\Users\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
  
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


load(file.path(models, "richnessmodel_correction.rds"))
load(file.path(models, "biomassmodel_full_correction.rds"))
load(file.path(models, "abundancemodel_full_correction.rds"))
# load(file.path(models, "fgrichnessmodel.rds"))



data_in <- "8_Data"



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

alldat$Study_Name <- as.character(alldat$Study_Name)
alldat$Study_Name[alldat$file == "4836_Hurisso2011"] <- "hurisso"
alldat$Study_Name <- as.factor(alldat$Study_Name)


## unique sites
length(unique(alldat$ID))
## 6931 # 9212
length(unique(alldat$Study_site))
## 6928 # 9207


length(unique(alldat$Study_Name)) # 228
length(unique(alldat$file)) # 180


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
countries <- sort(countries)

levels(alldat$country)


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

# 16
######################################################
## A UNIQUE DATASET
########################################################
uniqueDat <- alldat
alldat[,grep("scale", names(uniqueDat), ignore.case = TRUE)] <- NULL
uniqueDat <- uniqueDat[!duplicated(uniqueDat$Study_site), ]
## The scaled variables would be different for each dataframe, so removing them
nrow(uniqueDat)


summary(uniqueDat$PH)
nrow(uniqueDat) - 2504
summary(uniqueDat$CEC)
nrow(uniqueDat) -  6426


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


##########################################
# Number of sites with high number of months with snow
#########################################

table(uniqueDat$SnowMonths)
table(uniqueDat$SnowMonths_cat)


#####################################
# Max depth of weighted soil means
######################################
# I stupidly didn't code that sort of information
# So this will be a manual search (YAY!)
names(uniqueDat)[grep("_mean", names(uniqueDat))]

means1 <- as.vector(unique(uniqueDat$file[which(uniqueDat$PH_mean == "yes")]))
means2 <- as.vector(unique(uniqueDat$file[which(uniqueDat$CEC_mean == "yes")]))
means3 <- as.vector(unique(uniqueDat$file[which(uniqueDat$BaseSaturation_mean == "yes")]))
means4 <- as.vector(unique(uniqueDat$file[which(uniqueDat$OC_mean == "yes")]))
means5 <- as.vector(unique(uniqueDat$file[which(uniqueDat$SOM_mean == "yes")]))
means6 <- as.vector(unique(uniqueDat$file[which(uniqueDat$CN_mean == "yes")]))
means7 <- as.vector(unique(uniqueDat$file[which(uniqueDat$sand_silt_clay_mean == "yes")]))


means <- (c(means1, means2, means3, means4, means5, means6, means7))
unique(means) # 20

# [1] "000_KotanenUnpublished"   # 25cm       "000_Bescansa2010"    # 30cm (but not sure we weighted)          
# [3] "000_Virto2007"         # 1m         "3582_vanSchaik2014"            # 50cm
# [5] "000_FunDivEurope"     # 20cm          "7543_Bedano2016"      # 20cm         
# [7] "227_Richardson2015"    #10cm         "3331_Watmough2014"           # 10cm?  
# [9] "4130_Kuntz2013"         # 30cm        "1247_Rampazzo2001"            # 30cm 
# [11] "3225_Rozen1988"          #20cm       "384_Regulska2015"            30cm   
# [13] "5184_Li2010"    # 15cm (but not sure weighted)                "6736_Ammer2006"       30cm         
# [15] "4204_Uribe2012"         # 30cm        "4362_Castellanos-Navarrete2012"   #30 cm
# [17] "4414_Rahman2012"      #30cm          "4963_Ayuke2011"           # 30cm     
# [19] "5746_Ernst2009"         # 30cm        "4213_Virto2012"    # 30cm




## Where we don't have ph, are other soil variables coming from sampled data?
## 

noPh <- uniqueDat[is.na(uniqueDat$ph_new),]
summary(noPh$Organic_Carbon__percent)
summary(noPh$Silt__percent)
summary(noPh$Clay__percent)

haveClay <- droplevels(noPh[which(!(is.na(noPh$Silt__percent))),])
table(haveClay$Country)
