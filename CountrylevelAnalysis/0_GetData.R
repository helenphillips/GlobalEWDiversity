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

########################################################
# 5. Save the data
########################################################

write.csv(dat, file = file.path(data_out, "CountrySpeciesList.csv"), row.names = FALSE)


