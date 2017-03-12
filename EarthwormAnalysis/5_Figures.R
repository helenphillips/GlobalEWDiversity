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

library(lme4)
library(Hmisc)

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

sites$scalePH <-scale(sites$ph_new)
sites_habitat <- sites[sites$HabitatCover != "Unknown/Other",]
sites_lu <- droplevels(sites[sites$LandUse != "Unknown",])

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
### Interactions
plotInteraction(model = sp_habitat, Effect1 = "scalePH", Effect2 = "HabitatCover", 
                responseVar = "NumberofSpecies", seMultiplier = 1.96, 
                data = sites_habitat, cols = habitCols, legend.position = "topleft", 
                ylabel = "", xlabel = "")
  
plotInteraction(model = sp_lu, Effect1 = "scalePH", Effect2 = "LandUse", 
                responseVar = "NumberofSpecies", seMultiplier = 1.96, 
                data = sites_lu, cols = luCols, legend.position = "bottomleft", 
                ylabel = "", xlabel = "")



### Main effects
### Habitat cover
newdata <- with(sp_habitat@frame, expand.grid(scalePH = mean(scalePH, na.rm = TRUE), HabitatCover=levels(HabitatCover)))
newdata$NumberofSpecies <- predict(sp_habitat,newdata, re.form = NA)

mm <- model.matrix(terms(sp_habitat), newdata)

coef.names<-names(fixef(sp_habitat))
original.mm.names<-dimnames(mm)[[2]]
if(length(coef.names)!=length(original.mm.names)){mm<-mm[,dimnames(mm)[[2]]%in%coef.names]} 


pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(sp_habitat)), mm))
seMultiplier <- 1.96
newdata$upper <- newdata$NumberofSpecies + seMultiplier * sqrt(pvar1)
newdata$lower <- newdata$NumberofSpecies - seMultiplier * sqrt(pvar1)


newdata <- droplevels(newdata[newdata$HabitatCover != "Paddy field",])
newdata$cols <- ColourPicker(newdata$HabitatCover)
  newdata$cols2 <- paste("#", newdata$cols, sep="")

par(xpd=TRUE)
par(mar = c(3, 3, 2, 15))
plot(1:nrow(newdata) -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
      ylab = "log-Number of species", xlab = "Habitat cover", xaxt='n')
#xlim = c(min(newdata[,n],na.rm = TRUE), max(newdata[,n], na.rm = TRUE)), 

errbar(1:nrow(newdata), newdata$NumberofSpecies, newdata$upper, newdata$lower,
       add = TRUE, col = newdata$cols2, errbar.col = newdata$cols2, cex = 1.5)
legend(8.5,2,legend = newdata$HabitatCover, col = newdata$cols2, pch = 19, bty = "n", cex = 1.1)

# legend("bottomleft", legend = newdata$HabitatCover, col = newdata$cols2, lwd = 2, bty = "n")




### Land use
newdata <- with(sp_lu@frame, expand.grid(scalePH = mean(scalePH, na.rm = TRUE), LandUse=levels(LandUse)))
newdata$NumberofSpecies <- predict(sp_lu,newdata, re.form = NA)

mm <- model.matrix(terms(sp_lu), newdata)

coef.names<-names(fixef(sp_lu))
original.mm.names<-dimnames(mm)[[2]]
if(length(coef.names)!=length(original.mm.names)){mm<-mm[,dimnames(mm)[[2]]%in%coef.names]} 


pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(sp_lu)), mm))
seMultiplier <- 1.96
newdata$upper <- newdata$NumberofSpecies + seMultiplier * sqrt(pvar1)
newdata$lower <- newdata$NumberofSpecies - seMultiplier * sqrt(pvar1)

newdata$cols <- ColourPicker(newdata$LandUse)
newdata$cols2 <- paste("#", newdata$cols, sep="")

par(xpd=TRUE)
par(mar = c(3, 3, 2, 15))
plot(1:nrow(newdata) -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
     ylab = "log-Number of species", xlab = "Habitat cover", xaxt='n')
#xlim = c(min(newdata[,n],na.rm = TRUE), max(newdata[,n], na.rm = TRUE)), 

errbar(1:nrow(newdata), newdata$NumberofSpecies, newdata$upper, newdata$lower,
       add = TRUE, col = newdata$cols2, errbar.col = newdata$cols2, cex = 1.5)
legend(6.3,1.5,legend = newdata$LandUse, col = newdata$cols2, pch = 19, bty = "n", cex = 1.1)

# legend("bottomleft", legend = newdata$HabitatCover, col = newdata$cols2, lwd = 2, bty = "n")


