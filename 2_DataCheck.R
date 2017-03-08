
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
date <- "2017-03-08" ## TODO: Find a better way to do that

#################################################
# 4. Load in data
#################################################

sites <- read.csv(file.path(data_in, paste("sites_", date, ".csv", sep ="")))


#################################################
# 5. Get rid of studies with selected species 
#################################################

x <- tapply(sites$NumberofSpecies, sites$Study_Name, summary)

remove <- c("CoffeeFarms_Monoliths", "Malaysia_oilpalm", "Wu_Pheretima")

sites <- droplevels(sites[!(sites$Study_Name %in% remove),])

#################################################
# 6. Get rid of studies with selected species 
#################################################

write.csv(sites, file = file.path(data_out, paste("sites_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)


