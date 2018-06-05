########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

########################################################
# 2. Create folder if it doesn't exist to save data into
########################################################

if(!dir.exists("9_Data")){
  dir.create("9_Data")
}
data_out <- "9_Data"

data_in_fg <- "RevisedSpeciesNames"
data_in_spp <- "0_Data"
########################################################
# 3. Libraries
########################################################

library(googlesheets)
########################################################
# 4. Load data
########################################################

## Revised functional group datasets
fg <- read.csv(file.path(data_in_fg, "UniqueSpecies+FunctionalGroups_MJIB.csv"))
unwantedCol <- c("X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6", "X.7", "X.8", "X.9", "X.10", "X.11", "X.12")
fg <- fg[,!(names(fg) %in% unwantedCol)]

crucialCols <- c("original", "Revised", "Revised_fg")
fg <- fg[,(names(fg) %in% crucialCols)]

## Species level dataset
files <- list.files(file.path(data_in_spp))
files <- files[grep("species_", files)]


file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinspp <- loadin[grep("species_", loadin)]

spp <- read.csv(file.path(data_in_spp, loadinspp))

########################################################
# 5. Straight matches of binomial names
########################################################

spp <- merge(spp, fg, by.x = "SpeciesBinomial", by.y = "original", all.x = TRUE)
length(unique(spp$SpeciesBinomial)) # 306
length(unique(spp$Revised)) # 121


## There may be some where we have a binomial, but no match
spp2 <- spp[which(is.na(spp$Revised)),]
spp2 <- spp2[which(!(is.na(spp2$SpeciesBinomial))),]
spp2 <- as.data.frame(spp2$SpeciesBinomial)
names(spp2) <- "SpeciesBinomial"

## Some might match the revised version
spp2 <- unique(spp2)
spp2 <- merge(spp2, fg, by.x = "SpeciesBinomial", by.y = "Revised", all.x = TRUE)

complete <- spp2[which(!(is.na(spp2$Revised_fg))),]

for(species in 1:nrow(complete)){
  spp$Revised[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$SpeciesBinomial[species]
  spp$Revised_fg[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$Revised_fg[species]
}

## Back to just the missing ones
spp2 <- spp2[-which(complete.cases(spp2)),]
spp2$original <- NULL
spp2$Revised_fg <- NULL
## Some are definitely just slight spelling errors

for(fuzz in 1:nrow(spp2)){
  # print(fuzz)
  tri <- agrep(spp2$SpeciesBinomial[fuzz], fg$Revised, max.distance = 0.1)
  if(length(tri) > 0){
    tri <- tri[1]
    spp2$FuzzyName[fuzz] <- as.character(fg$Revised[tri])
    spp2$FuzzyFg[fuzz] <- as.character(fg$Revised_fg[tri])
  }
  if(length(tri) == 0){
    spp2$FuzzyName[fuzz] <- NA
    spp2$FuzzyFg[fuzz] <- NA
  }
}

complete <- spp2[which(!(is.na(spp2$FuzzyFg))),]
## The first line is not a good match
complete <- complete[-1,]

for(species in 1:nrow(complete)){
  spp$Revised[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$FuzzyName[species]
  spp$Revised_fg[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$FuzzyFg[species]
}

### Some may be in drilobase

#################################################
# 9. DriloBASE data
#################################################

drilobase <- "DriloBASE_uniqueSpecies"
drilobase <- gs_title(drilobase)
drilo <- as.data.frame(gs_read(drilobase, ws = "Sheet1"))


## Back to just the missing ones
spp2[4,c(2:3)] <- NA

spp2 <- spp2[-which(complete.cases(spp2)),]
spp2$FuzzyName <- NULL
spp2$FuzzyFg <- NULL

for(fuzz in 1:nrow(spp2)){
  # print(fuzz)
  tri <- agrep(spp2$SpeciesBinomial[fuzz], drilo$name, max.distance = 0.1)
  if(length(tri) > 0){
    tri <- tri[1]
    spp2$FuzzyName[fuzz] <- as.character(drilo$name[tri])
    spp2$FuzzyFg[fuzz] <- as.character(drilo$ecologicalCategory[tri])
  }
  if(length(tri) == 0){
    spp2$FuzzyName[fuzz] <- NA
    spp2$FuzzyFg[fuzz] <- NA
  }
}

## row 4 and 124 are not good matches
spp2[4,c(2:3)] <- NA
spp2[which(spp2$SpeciesBinomial == "Prosellodrilus praticola"), c(2:3)] <- NA

complete <- spp2[which(!(is.na(spp2$FuzzyName))),]

spp$Revised <- as.character(spp$Revised)
spp$Revised_fg <- as.character(spp$Revised_fg)

for(species in 1:nrow(complete)){
  spp$Revised[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$FuzzyName[species]
  spp$Revised_fg[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$FuzzyFg[species]
}

## Back to just the missing ones

spp2 <- spp2[which(is.na(spp2$FuzzyName)),] # 49 rows
spp2$FuzzyName <- NULL
spp2$FuzzyFg <- NULL


## 49 left
spp2$HelenMatch <- NA
spp2$HelenFG <- NA

spp2[spp2$SpeciesBinomial == "Allobophora caliginosa",2:3] <- c("Aporrectodea caliginosa", "Endogeic")
spp2[spp2$SpeciesBinomial == "Aporrectodea caliginosa caliginosa",2:3] <- c("Aporrectodea caliginosa caliginosa", "Endogeic")
spp2[spp2$SpeciesBinomial == "Aporrectodea oliveirae oliveirae",2:3] <- c("Aporrectodea oliveirae oliveirae", "NA")
spp2[spp2$SpeciesBinomial == "Aporrectodea oliveirae trigoae",2:3] <- c("Aporrectodea oliveirae trigoae", "NA")
spp2[spp2$SpeciesBinomial == "Aporrectodea rosea bimastoides",2:3] <- c("Aporrectodea rosea bimastoides", "Endogeic")
spp2[spp2$SpeciesBinomial == "Aporrectodea rosea rosea",2:3] <- c("Aporrectodea rosea rosea", "Endogeic")
spp2[spp2$SpeciesBinomial == "Bimastus tenius",2:3] <- c("Dendrodrilus rubidus tenuis", "Epigeic")
spp2[spp2$SpeciesBinomial == "Dendrobaena mammalis",2:3] <- c("Satchellius mammalis", "Epigeic")
spp2[spp2$SpeciesBinomial == "Dendrobaena pygmea cognettii",2:3] <- c("Dendrobaena pygmea cognettii", "Epigeic")
spp2[spp2$SpeciesBinomial == "Dendrobaena subrubicunda",2:3] <- c("Dendrodrilus rubidus subrubicundus", "Epigeic")
spp2[spp2$SpeciesBinomial == "Dendrobaena veneta",2:3] <- c("Eisenia veneta", "Epigeic")
spp2[spp2$SpeciesBinomial == "Eminoscolex violaceus",2:3] <- c("Eminoscolex violaceus", "NA")
spp2[spp2$SpeciesBinomial == "Eudrilus buettneri ifensis",2:3] <- c("Eudrilus pallidus", "NA")
spp2[spp2$SpeciesBinomial == "Heraclescolex moebii",2:3] <- c("Aporrectodea moebii", "NA")
spp2[spp2$SpeciesBinomial == "Iberoscolex albolineatus",2:3] <- c("Iberoscolex albolineatus", "NA")
spp2[spp2$SpeciesBinomial == "Karpatodinariona altimontana",2:3] <- c("Allolobophora altimontana", "NA")
spp2[spp2$SpeciesBinomial == "Kritodrilus pseudorroseus",2:3] <- c("Iberoscolex pseudorrosea", "Endogeic")
spp2[spp2$SpeciesBinomial == "Metaphire inflata cai",2:3] <- c("Metaphire cai", "NA")
spp2[spp2$SpeciesBinomial == "Microeophila nematogena",2:3] <- c("Perelia nematogena", "Endogeic")
spp2[spp2$SpeciesBinomial == "Octolasion tyrtaeum tyrtaeum",2:3] <- c("Octolasion tyrtaeum tyrtaeum", "Endogeic")
spp2[spp2$SpeciesBinomial == "Prosellodrilus amplisetosus amplisetosus",2:3] <- c("Prosellodrilus amplisetosus amplisetosus", "NA")
spp2[spp2$SpeciesBinomial == "Prosellodrilus praticola",2:3] <- c("Prosellodrilus praticola", "NA")
spp2[spp2$SpeciesBinomial == "Righiodrilus gurupi",2:3] <- c("Righiodrilus gurupi", "NA")

