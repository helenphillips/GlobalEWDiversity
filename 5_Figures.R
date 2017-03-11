########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
source(file.path("Functions", "Plots.R"))
source(file.path("Functions", "ColourPicker.R"))

#################################################
# 2. Loading in variables
#################################################

data_in <-"3_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates)
loadin <- files[grep(date, files)]

rm(files)
rm(date)

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

#################################################
# 3. Load in data
#################################################

sites <- read.csv(file.path(data_in, loadin))
rm(loadin)


#################################################
# 4. Load in models
#################################################

models <- "Models"

load(file.path(models, "sp_habitat.rds"))
load(file.path(models, "sp_landuse.rds"))

#################################################
# 5. Pick colors for figures
#################################################
habitCols <- ColourPicker(sites_habitat$HabitatCover)
luCols <- ColourPicker(sites_lu$LandUse)
#################################################
# 6. Figures
#################################################

plotInteraction(model = sp_habitat, Effect1 = "scalePH", Effect2 = "HabitatCover", 
                responseVar = "NumberofSpecies", seMultiplier = 1.96, 
                data = sites_habitat, cols = habitCols, legend.position = "topleft", 
                ylabel = "", xlabel = "")
  
plotInteraction(model = sp_lu, Effect1 = "scalePH", Effect2 = "LandUse", 
                responseVar = "NumberofSpecies", seMultiplier = 1.96, 
                data = sites_lu, cols = luCols, legend.position = "bottomleft", 
                ylabel = "", xlabel = "")

