
########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

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

data_in <-"0_Data"
files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
loadin <- files[grep(date, files)]
loadin <- loadin[grep("sites_", loadin)]


#################################################
# 4. Load in data
#################################################

sites <- read.csv(file.path(data_in, loadin))


#################################################
# 5. Get rid of studies with selected species 
#################################################

x <- tapply(sites$NumberofSpecies, sites$Study_Name, summary)

remove <- c("CoffeeFarms_Monoliths", "Malaysia_oilpalm", "Wu_Pheretima")

sites <- droplevels(sites[!(sites$Study_Name %in% remove),])


#################################################
# MAke some factor levels consistent
#################################################

levels(sites$HabitatCover)[levels(sites$HabitatCover) == "Unknown"] <- "Unknown/Other"


#################################################
# 6. Save file
#################################################

write.csv(sites, file = file.path(data_out, paste("sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)


