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
data_in_spp <- "1_Data"
########################################################
# 3. Libraries
########################################################

library(googlesheets)
########################################################
# 4. Load data
########################################################
## Species level dataset
files <- list.files(file.path(data_in_spp))
files <- files[grep("UniqueSpecies_", files)]

file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinspp <- loadin[grep("UniqueSpecies_", loadin)]

spp <- read.csv(file.path(data_in_spp, loadinspp))


## Already Revised functional group datasets
fg <- read.csv(file.path(data_in_fg, "UniqueSpecies+FunctionalGroups_MJIB.csv"))
unwantedCol <- c("X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6", "X.7", "X.8", "X.9", "X.10", "X.11", "X.12")
fg <- fg[,!(names(fg) %in% unwantedCol)]

########################################################
# 4.5 Edit data
########################################################
# Update unique species with already revised
crucialCols <- c("original", "Revised", "Revised_fg", "Revised_Authority", "Notes")
fg <- fg[,(names(fg) %in% crucialCols)]

spp <- merge(spp, fg, by.x = "SpeciesBinomial", by.y = "original", all.x = TRUE)
 
########################################################
# 5. Straight matches of binomial names
########################################################

length(unique(spp$SpeciesBinomial)) # 304
length(unique(spp$Revised)) # 121

## There may be some where we have a binomial, but no match
spp2 <- spp[which(is.na(spp$Revised)),]
spp2[,which(names(spp2) %in% c("Revised", "Revised_fg", "Revised_Authority", "Notes"))] <- NULL
## But might be because they are already correct
spp2 <- merge(spp2, fg, by.x = "SpeciesBinomial", by.y = "Revised", all.x = TRUE)

complete <- spp2[which(!(is.na(spp2$Revised_Authority))),]

## Add to SppList


for(species in 1:nrow(complete)){
  spp$Revised[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$SpeciesBinomial[species]
  spp$Revised_fg[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$Revised_fg[species]
  spp$Revised_Authority[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$Revised_Authority[species]
  spp$Notes[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$Notes[species]
}

## Back to just the missing ones
spp2 <- spp2[-which(complete.cases(spp2$Revised_Authority)),]
spp2[,which(names(spp2) %in% c("Revised", "Revised_fg", "Revised_Authority", "Notes"))] <- NULL

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

spp$Revised_Authority <- as.character(spp$Revised_Authority)
spp$Notes <- as.character( spp$Notes)

for(species in 1:nrow(complete)){
  spp$Revised[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$FuzzyName[species]
  spp$Revised_fg[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$FuzzyFg[species]
  spp$Revised_Authority[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- "Helen - spelling errors"
  spp$Notes[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- "Used species names already corrected by EW experts"
}



### Some may be in drilobase
#################################################
# 8. Maria had already checked species names in some of the datasets
#################################################
## Back to just the missing ones
spp2[4,6:8] <- NA

spp2 <- spp2[-which(complete.cases(spp2$FuzzyName)),]
spp2$original <- NULL
spp2$FuzzyName <- NULL
spp2$FuzzyFg <- NULL

inde <- grep("000_Briones1991", spp2$file)
inde <- c(inde, grep("000_Trigo1987", spp2$file))
inde <- c(inde, grep("7601_GutierrezLopez2016", spp2$file))
inde <- unique(inde)
# 5 13 14 16 17 27 46 58 62 63 65 93 94 97 98 12 15 18 26 29 31 33 37 

complete <- spp2[inde,]

spp$Revised <- as.character(spp$Revised)
spp$Revised_fg <- as.character(spp$Revised_fg)


for(species in 1:nrow(complete)){
  spp$Revised[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- as.character(complete$SpeciesBinomial[species])
  spp$Revised_fg[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- as.character(complete$FunctionalGroup[species])
  spp$Revised_Authority[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- "Briones"
  spp$Notes[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- "Checked prior to data entry"
}


#################################################
# 9. DriloBASE data
#################################################

drilobase <- "DriloBASE_uniqueSpecies"
drilobase <- gs_title(drilobase)
drilo <- as.data.frame(gs_read(drilobase, ws = "Sheet1"))


## Back to just the missing ones

spp2 <- spp2[-inde,]


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

spp2[,c(1, 6)]
## row 4 and 124 are not good matches
spp2[which(spp2$SpeciesBinomial == "Allolobophora parva"), 6:7] <- NA
spp2[which(spp2$SpeciesBinomial == "Prosellodrilus praticola"), 6:7] <- NA

complete <- spp2[which(!(is.na(spp2$FuzzyName))),]


for(species in 1:nrow(complete)){
  spp$Revised[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$FuzzyName[species]
  spp$Revised_fg[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$FuzzyFg[species]
  spp$Revised_Authority[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- "Helen - Drilobase"
  spp$Notes[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- "Fuzzy match to Drilobase"
  
}

## Back to just the missing ones

spp2 <- spp2[which(is.na(spp2$FuzzyName)),] # 49 rows
spp2$FuzzyName <- NULL
spp2$FuzzyFg <- NULL


## 43 left
spp2$HelenMatch <- NA
spp2$HelenFG <- NA

spp2[spp2$SpeciesBinomial == "Allobophora caliginosa",6:7] <- c("Aporrectodea caliginosa", "Endogeic")
spp2[spp2$SpeciesBinomial == "Aporrectodea caliginosa caliginosa",6:7] <- c("Aporrectodea caliginosa caliginosa", "Endogeic")
# spp2[spp2$SpeciesBinomial == "Aporrectodea oliveirae oliveirae",6:7] <- c("Aporrectodea oliveirae oliveirae", "NA")
# spp2[spp2$SpeciesBinomial == "Aporrectodea oliveirae trigoae",6:7] <- c("Aporrectodea oliveirae trigoae", "NA")
spp2[spp2$SpeciesBinomial == "Aporrectodea rosea bimastoides",6:7] <- c("Aporrectodea rosea bimastoides", "Endogeic")
spp2[spp2$SpeciesBinomial == "Aporrectodea rosea rosea",6:7] <- c("Aporrectodea rosea rosea", "Endogeic")
spp2[spp2$SpeciesBinomial == "Bimastus tenius",6:7] <- c("Dendrodrilus rubidus tenuis", "Epigeic")
spp2[spp2$SpeciesBinomial == "Dendrobaena mammalis",6:7] <- c("Satchellius mammalis", "Epigeic")
spp2[spp2$SpeciesBinomial == "Dendrobaena pygmea cognettii",6:7] <- c("Dendrobaena pygmea cognettii", "Epigeic")
spp2[spp2$SpeciesBinomial == "Dendrobaena subrubicunda",6:7] <- c("Dendrodrilus rubidus subrubicundus", "Epigeic")
# spp2[spp2$SpeciesBinomial == "Dendrobaena veneta",6:7] <- c("Eisenia veneta", "Epigeic")
spp2[spp2$SpeciesBinomial == "Eminoscolex violaceus",6:7] <- c("Eminoscolex violaceus", "NA")
spp2[spp2$SpeciesBinomial == "Eudrilus buettneri ifensis",6:7] <- c("Eudrilus pallidus", "NA")
spp2[spp2$SpeciesBinomial == "Heraclescolex moebii",6:7] <- c("Aporrectodea moebii", "NA")
spp2[spp2$SpeciesBinomial == "Iberoscolex albolineatus",6:7] <- c("Iberoscolex albolineatus", "NA")
spp2[spp2$SpeciesBinomial == "Karpatodinariona altimontana",6:7] <- c("Allolobophora altimontana", "NA")
spp2[spp2$SpeciesBinomial == "Kritodrilus pseudorroseus",6:7] <- c("Iberoscolex pseudorrosea", "Endogeic")
spp2[spp2$SpeciesBinomial == "Metaphire inflata cai",6:7] <- c("Metaphire cai", "NA")
spp2[spp2$SpeciesBinomial == "Microeophila nematogena",6:7] <- c("Perelia nematogena", "Endogeic")
spp2[spp2$SpeciesBinomial == "Octolasion tyrtaeum tyrtaeum",6:7] <- c("Octolasion tyrtaeum tyrtaeum", "Endogeic")
spp2[spp2$SpeciesBinomial == "Prosellodrilus amplisetosus amplisetosus",6:7] <- c("Prosellodrilus amplisetosus amplisetosus", "NA")
spp2[spp2$SpeciesBinomial == "Prosellodrilus praticola",6:7] <- c("Prosellodrilus praticola", "NA")
spp2[spp2$SpeciesBinomial == "Righiodrilus gurupi",6:7] <- c("Righiodrilus gurupi", "NA")

complete <- spp2[which(!(is.na(spp2$HelenMatch))),]

for(species in 1:nrow(complete)){
  spp$Revised[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$HelenMatch[species]
  spp$Revised_fg[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- complete$HelenFG[species]
  spp$Revised_Authority[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- "Helen - Drilobase Synonyms"
  spp$Notes[which(spp$SpeciesBinomial == complete$SpeciesBinomial[species])] <- "Manual match to Drilobase"
  
}

### Back to just missing ones
spp2 <- spp2[which(is.na(spp2$HelenMatch)),] # 25 rows
spp2$HelenMatch <- NULL
spp2$HelenFG <- NULL


## But the spp file should now be everything filled in
#################################################
# 11. Save Data
#################################################

write.csv(spp, file = file.path(data_out, paste("Unique_Species_toSend", Sys.Date(), ".csv", sep = "")), row.names = FALSE)
