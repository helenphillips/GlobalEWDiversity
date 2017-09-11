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
# habitCols <- ColourPicker(sites_habitat$HabitatCover)
# luCols <- ColourPicker(sites_lu$LandUse)
luMgmtCols <- ColourPicker(abundance$LU_Mgmt)
#################################################
# 6. Figures
#################################################
### Interactions
plotInteraction(model = a1, Effect1 = "scalePH", Effect2 = "LU_Mgmt", 
                responseVar = "logAbundance", seMultiplier = 1.96, 
                data = abundance, cols = luMgmtCols, legend.position = "topleft", 
                ylabel = "", xlabel = "")
  
#plotInteraction(model = sp_lu, Effect1 = "scalePH", Effect2 = "LandUse", 
#                responseVar = "NumberofSpecies", seMultiplier = 1.96, 
#                data = sites_lu, cols = luCols, legend.position = "bottomleft", 
#                ylabel = "", xlabel = "")



### Main effects
newdata <- with(a1@frame, expand.grid(scalePH = mean(scalePH, na.rm = TRUE), LU_Mgmt=levels(LU_Mgmt)))
newdata$logAbundance <- predict(a1,newdata, re.form = NA)

mm <- model.matrix(terms(a1), newdata)

coef.names<-names(fixef(a1))
original.mm.names<-dimnames(mm)[[2]]
if(length(coef.names)!=length(original.mm.names)){mm<-mm[,dimnames(mm)[[2]]%in%coef.names]} 


pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(a1)), mm))
seMultiplier <- 1.96
newdata$upper <- newdata$logAbundance + seMultiplier * sqrt(pvar1)
newdata$lower <- newdata$logAbundance - seMultiplier * sqrt(pvar1)

#newdata$cols <- ColourPicker(newdata$HabitatCover)
#newdata$cols2 <- paste("#", newdata$cols, sep="")

# pdf(file = file.path(figures, "HabitatCover.pdf"), width = 11)
#par(xpd=TRUE)
#par(mar = c(2, 4.5, 2, 16))
plot(1:nrow(newdata) -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
     ylab = "log-Abundance", xlab = "LU_Mgmt", xaxt='n', cex.lab = 1.5)
#xlim = c(min(newdata[,n],na.rm = TRUE), max(newdata[,n], na.rm = TRUE)), 

errbar(1:nrow(newdata), newdata$logAbundance, newdata$upper, newdata$lower,
       add = TRUE, col = luMgmtCols, errbar.col = luMgmtCols, cex = 1.5)
text(1:8, 2.3, paste("n =", table(a1@frame$LU_Mgmt)))
# dev.off()
