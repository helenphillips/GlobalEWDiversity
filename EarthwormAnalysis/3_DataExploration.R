########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



#################################################
# 1. Loading libraries
#################################################
library(maptools)
library(maps)
#################################################
# 2. Loading in variables
#################################################

data_in <-"2_Data"
date <- "2017-03-08" ## TODO: Find a better way to do that

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"
#################################################
# 3. Load in data
#################################################

sites <- read.csv(file.path(data_in, paste("sites_", date, ".csv", sep ="")))


#################################################
# 4. Basic stats
#################################################

length(unique(sites$file)) ## 13 papers
length(unique(sites$Study_Name)) ## 15 studies

length(unique(sites$Country))## 11 Countries

#################################################
# 5. Create Map
#################################################

coord<-aggregate(cbind(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees), list(sites$Study_Name), mean)
## Three don't have coordinates yet
coord <- coord[complete.cases(coord),]


nocoords <- c("FalcoWorms", "Henneron", "Regulska_poland")
unique(sites$Country[sites$Study_Name %in% nocoords])
# Argentina France    Poland
missing <- data.frame(country = c("Argentina", "France",  "Poland"), Long = c(-63.6167, 2.2137, 19.1451), 
                      Lat = c(-38.4161, 46.2276, 51.9194))
  


coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")

missingSPDF<-SpatialPointsDataFrame(missing[,2:3], data.frame(missing[,1:3]))
proj4string(missingSPDF)<-CRS("+proj=longlat")


pdf(file = file.path(figures, "Map.pdf"), height = 4)
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)
points(missingSPDF, col="red", bg="red", cex= 1, pch=19)
dev.off()



######################################################
## 
######################################################


tapply(sites$PH, sites$Study_Name, summary)
## 4 with none 

tapply(sites$Organic_Carbon__percent, sites$Study_Name, summary)
## 11 with none

tapply(sites$Soil_Organic_Matter__percent, sites$Study_Name, summary)
## Only three with

tapply(sites$C.N_ratio, sites$Study_Name, summary)
## Only three with

table(sites$LandUse, sites$Study_Name)
## All land uses have comparisons

table(sites$HabitatCover, sites$Study_Name)
