########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

########################################################
# 2. Create folder if it doesn't exist to save data into
########################################################

if(!dir.exists("1_Data")){
  dir.create("1_Data")
}
data_out <- "1_Data"


########################################################
# 3. Libraries
########################################################

source("Functions/FormatData.R")
library(dplyr)
library(googlesheets)
x <- gs_ls() ## Authentication

#################################################
# 4. Loading in variables
#################################################

data_in <-"0_Data"
files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
loadin <- files[grep(date, files)]
loadin <- loadin[grep("species_", loadin)]


#################################################
# 5. Load in data
#################################################

dat <- read.csv(file.path(data_in, loadin))

#################################################
# 6. Quick investigation
#################################################

length(unique(dat$SpeciesBinomial)) ## 180

table(dat$Functional_Type) ## Only unknown for 310

## Check that all are binomials
which(
  (sapply(gregexpr("\\W+", dat$SpeciesBinomial), length) + 1) 
  < 2) ## This counts the number of words in teh string
#################################################
# 7. Create dataframe for new morphospecies
#################################################

result <- dat %>%
   group_by(SpeciesBinomial) %>%
   summarise(strings = toString(unique(Functional_Type)))

spp <- as.data.frame(result)
spp$strings <- gsub("NA, ", "", spp$strings)
spp$strings <- gsub(", NA", "", spp$strings)
spp$strings <- gsub("Unknown, ", "", spp$strings)
spp$strings <- gsub(", Unknown", "", spp$strings)

names(spp)[2] <- "FunctionalGroup"

## Remove line that isn't a species
spp <- spp[-which(is.na(spp$SpeciesBinomial)),]

## Remove Microdriles and Megadriles
spp <- spp[-grep("microdrile", spp$SpeciesBinomial, ignore.case = TRUE),]
spp <- spp[-grep("Megadrile", spp$SpeciesBinomial, ignore.case = TRUE)]


#################################################
# 9. DriloBASE data
#################################################

drilobase <- "DriloBASE_uniqueSpecies"
drilobase <- gs_title(drilobase)
drilo <- as.data.frame(gs_read(drilobase, ws = "Sheet1"))

###################################################
# 10. Map our dataset
###############################################

s2 <- adist(spp$SpeciesBinomial, drilo$name)
i2 <- apply(s2, 1, which.min) ## Gives the index
spp$Drilobase <- drilo$name[i2]

cols <-c("SpeciesBinomial", "FunctionalGroup", "Drilobase", "Author of species", "ecologicalCategory")
spp <- (merge(spp, drilo, by.x = "Drilobase", by.y = "name", all.x = TRUE))
spp <- spp [,names(spp) %in% cols]
spp <- spp[,c(2, 3, 1, 4, 5)]
names(spp) <- c("original", "original_fg", "drilobase", "Authority of species", "drilobase_fg")

###################################################
# 111. Add blank columns for people to fill in
###############################################


spp$Revised <- NA
spp$Revised_fg <- NA
spp$Revised_Authority <- NA
spp$sWormMember <- NA

###################################################
# 11. Save the data
###############################################

write.csv(spp, file.path(data_out, paste("UniqueSpecies_", Sys.Date(), ".csv", sep ="")), row.names = FALSE)

#Also need to save to google sheets, which is where everyone will edit
## This has been done, and until you download what people have filled in
## DO NOT run this again
## Also, it may not work

# output <- "UniqueSpecies+FunctionalGroups"
# output <- gs_title(output)
# output <- (gs_gs(output))
# 
# t <- try(output <- output %>% 
#   gs_ws_new(ws_title = "UniqueSpecies+FunctionalGroups", input = spp,
#             trim = TRUE, verbose = FALSE), silent = TRUE)
# 
# if(class(t) == "try-error"){ ## Googlesheet already exists, need to delete old one
#   output <- output %>% 
#     gs_ws_delete(ws = "UniqueSpecies+FunctionalGroups")
#   
#   ## And upload new one
#   output <- output %>% 
#     gs_ws_new(ws_title = "UniqueSpecies+FunctionalGroups", input = spp,
#               trim = TRUE, verbose = FALSE)
#   
# }
# 
# output %>% 
#   gs_read(ws = 1)

