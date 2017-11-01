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

length(unique(dat$SpeciesBinomial)) ## 187

sort(levels(dat$SpeciesBinomial))


#################################################
# 7. Create dataframe for new morphospecies
#################################################

spp <- data.frame(original = sort(levels(dat$SpeciesBinomial)))

#################################################
# 8. Tally of species names 
#################################################
sppN <- as.data.frame(table(dat$SpeciesBinomial))

spp <- (merge(spp, sppN, by.x = "original", by.y = "Var1"))

# Take first hundred most predominant names
popular <- spp$original[order(spp$Freq, decreasing = TRUE)]
popular <- popular[1:30]

s <- adist(spp$original, popular)
i <- apply(s, 1, which.min) ## Gives the index
spp$Proposed <- popular[i]

#################################################
# 9. DriloBASE data
#################################################

drilobase <- "DriloBASE_uniqueSpecies"
drilobase <- gs_title(drilobase)
drilo <- as.data.frame(gs_read(drilobase, ws = "Sheet1"))

###################################################
# 10. Map our dataset
###############################################

s2 <- adist(spp$original, drilo$name)
i2 <- apply(s2, 1, which.min) ## Gives the index
spp$Drilobase <- drilo$name[i2]

spp <- (merge(spp, drilo, by.x = "Drilobase", by.y = "name", all.x = TRUE))[,c(2, 3, 4, 1, 10)]
names(spp)[5] <- "drilobase_fg"


###################################################
# 11. Save the data
###############################################

write.csv(spp, file.path(data_out, "UniqueSpecies.csv"), row.names = FALSE)
