########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\CountryLevelAnalysis\\")
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

require(RJSONIO)    

########################################################
# 4. Load data and sort out data
########################################################
## Species by country
json_file <- "C://Users//hp39wasi//Google Drive//sWorm//CountryLevelAnalysis//dataDriloBASE.json"

json_file <- fromJSON(json_file)

dat <- data.frame(Country = NA, species = NA)

for(i in 1:length(json_file)){
  t <- json_file[[i]]
  name <- t[1]
  species <- t[-1]
  
  df <- data.frame(Country = rep(NA, length(species)), species = rep(NA, length(species)))
  df[,1] <- rep(name, length(species))
  df[,2] <- species
  
  dat <- rbind(dat, df)
  
}

## Remove first row
dat <- dat[-1,]


### Species information
json_file <- "C://Users//hp39wasi//Google Drive//sWorm//CountryLevelAnalysis//speciesData.json"
json_file <- fromJSON(json_file)

species <- data.frame(name = NA,family = NA, genus = NA, geographicalOrigin = NA, distributionStatus = NA, ecologicalCategory = NA)

for(i in 1:length(json_file)){
  t <- json_file[[i]]
  species[i,'name'] <- t['name']
  species[i,'family'] <- t['family']
  species[i,'genus'] <- t['genus']
  species[i,'geographicalOrigin'] <- t['geographicalOrigin']
  species[i,'distributionStatus'] <- t['distributionStatus']
  species[i,'ecologicalCategory'] <- t['ecologicalCategory']
  
}

species$name <- sub("_", " ", species$name)
########################################################
# 5. Save the data
########################################################

write.csv(dat, file = file.path(data_out, "CountrySpeciesList.csv"), row.names = FALSE)
write.csv(species, file = file.path(data_out, "SpeciesList.csv"), row.names = FALSE)



########################################################
# 6.Amalgamate data
########################################################
c <- read.csv(file.path(data_out, "CountrySpeciesList.csv"))
s <- read.csv(file.path(data_out, "SpeciesList.csv"))

CountryAndSpecies <- merge(c, s, by.x = "species", by.y = "name", all.x = TRUE)


CountryAndSpecies <- CountryAndSpecies[with(CountryAndSpecies, order(Country, species)), ]
write.csv(CountryAndSpecies, file = file.path(data_out, "CountryAndSpeciesList.csv"), row.names = FALSE)





### Just want the list of countries for Alberto

countries <- data.frame(Country = unique(c$Country))
write.csv(countries, file = file.path(data_out, "CountryList.csv"), row.names = FALSE)
