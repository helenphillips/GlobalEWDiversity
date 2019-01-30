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
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinsites <- loadin[grep("species_", loadin)]

loadinbib <- loadin[grep("Metadata_", loadin)]
#################################################
# 5. Load in data
#################################################

dat <- read.csv(file.path(data_in, loadinsites))
bib <- read.csv(file.path(data_in, loadinbib))

#################################################
# 6. Quick investigation
#################################################

length(unique(dat$SpeciesBinomial)) ##  306

table(dat$Functional_Type) ## Only unknown for  821

## Check that all are binomials
which(
  (sapply(gregexpr("\\W+", dat$SpeciesBinomial), length) + 1) 
  < 2) ## This counts the number of words in teh string
#################################################
# 6. Merge bib info
#################################################

bib <- bib[,which(names(bib) == "file" | names(bib) == "DataProvider_Surname")]

dat <- merge(dat, bib, by.x = "file.x", by.y = "file")

#################################################
# 7. Create dataframe for new morphospecies
#################################################

result <- dat %>%
   group_by(SpeciesBinomial) %>%
   summarise(fg = toString(unique(Functional_Type)), 
             Country = toString(unique(Country)), 
             file = toString(unique(file.x)),
             provider = toString(unique(DataProvider_Surname)))

spp <- as.data.frame(result)
## Remove line that isn't a species
spp <- spp[-which(is.na(spp$SpeciesBinomial)),]
spp$fg <- gsub("NA, ", "", spp$fg)
spp$fg <- gsub(", NA", "", spp$fg)
spp$fg <- gsub("Unknown, ", "", spp$fg)
spp$fg <- gsub(", Unknown", "", spp$fg)

names(spp)[2] <- "FunctionalGroup"



## Remove Microdriles and Megadriles
spp <- spp[-grep("microdrile", spp$SpeciesBinomial, ignore.case = TRUE),]
spp <- spp[-grep("Megadrile", spp$SpeciesBinomial, ignore.case = TRUE)]


#################################################
# 9. DriloBASE data
#################################################
# 
# drilobase <- "DriloBASE_uniqueSpecies"
# drilobase <- gs_title(drilobase)
# drilo <- as.data.frame(gs_read(drilobase, ws = "Sheet1"))
# 
# ###################################################
# # 10. Map our dataset
# ###############################################
# 
# s2 <- adist(spp$SpeciesBinomial, drilo$name)
# i2 <- apply(s2, 1, which.min) ## Gives the index
# spp$Drilobase <- drilo$name[i2]
# 
# cols <-c("SpeciesBinomial", "FunctionalGroup", "Country", "file", "provider", 
#          "Drilobase", "Author of species", "ecologicalCategory")
# spp <- (merge(spp, drilo, by.x = "Drilobase", by.y = "name", all.x = TRUE))
# spp <- spp [,names(spp) %in% cols]
# spp <- spp[,c(2, 3, 4, 5, 6, 1, 7, 8)]
# names(spp) <- c("original", "original_fg", "Country", "PaperID", "dataProvider", "drilobase", "Authority of species", "drilobase_fg")
# 
# ###################################################
# # 111. Add blank columns for people to fill in
# ###############################################
# 
# 
# spp$Revised <- NA
# spp$Revised_fg <- NA
# spp$Revised_Authority <- NA
# spp$sWormMember <- NA

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

######################################
## Match to any previous Google sheet
######################################

output <- "UniqueSpecies+FunctionalGroups"
output <- gs_title(output)
current <- as.data.frame(gs_read(output, ws = "sheet1"))

## names(current) <- gsub("\\.y", "", names(current))
## names(current) <- gsub("\\.x", "", names(current))


spp <- merge(spp, current, by = "original")
keep <- c("original","original_fg.x","Country","PaperID", "dataProvider.x", 
          "drilobase.x", "Authority of species.x", "drilobase_fg.x", "Revised.y",
          "Revised_fg.y","Revised_Authority.y","sWormMember.y", "X12") 
spp <- spp[,which(names(spp) %in% keep)]

spp$Revised.y[which(spp$sWormMember.y %in% c("Patrick Lavelle", "Patrick") & is.na(spp$Revised.y))] <- as.character(spp$original[which(spp$sWormMember.y %in% c("Patrick Lavelle", "Patrick") & is.na(spp$Revised.y))])

names(spp) <- gsub("\\.y", "", names(spp))
names(spp) <- gsub("\\.x", "", names(spp))

write.csv(spp, file.path(data_out, paste("UniqueSpecies-Revised_", Sys.Date(), ".csv", sep ="")), row.names = FALSE)
