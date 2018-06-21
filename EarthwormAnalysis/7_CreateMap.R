if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}


if(!dir.exists("Figures")){
  dir.create("Figures")
}
figures <- "Figures"

## Create map
library(raster)
library(RColorBrewer)


# spR_finalraster <- raster(file.path(soil_GLs, "spRFinalRaster.tif"))
spR_finalraster <- raster("I:\\sDiv-postdocs-work\\Phillips\\sWorm\\SpatialAnalysis\\ProcessedLayers\\April2018\\spRFinalRasterminusESA.tif")



# summary(spR_finalraster)
# 0.26, 2.19  (min, max)

blubrks <- seq( -0.05,1.175, length.out = 100)
blues <- rev(brewer.pal(9, "Blues"))
b.cols <- colorRampPalette(blues, space = "rgb")


redbrks <-  seq(1.176, 2.33, length.out = 100)
reds <- brewer.pal(9, "Reds")
r.cols <- colorRampPalette(reds, space = "rgb")

png(file.path(figures, "SpeciesRichness.png"),width=17.5,height=8.75,units="cm",res=1200)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))
image(spR_finalraster, col=b.cols(length(blubrks)-1), breaks=blubrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(spR_finalraster, col=r.cols(length(redbrks)-1), add = TRUE, breaks=redbrks, xaxt="n", yaxt="n", ylab="", xlab="")
#hist(1:10, axes = FALSE)
par(mar=c(1,13,1,13))
blu_scale <- rep(b.cols(length(blubrks)-1), each = 2) # 200 total
red_scale <-rep(r.cols(length(redbrks)-1), each = 2)



barplot(rep(1, 396), col = c(blu_scale, red_scale), border =c(blu_scale, red_scale), axes = FALSE )
mtext("1", at = 20, cex = 1)
# mtext("6", at = 350, cex = 0.5)
mtext("8", at = 455, cex = 1)

dev.off()

rm(sr)
gc()


####### Biomass

# regions <- c("africa", "asia", "europe", "latin_america", "north_america", "west_asia")

results <- "I:\\sDiv-PostDocs-Work\\Phillips\\sWorm\\SpatialAnalysis\\Results\\Biomass"

 regions <- c("asia", "europe", "latin_america", "north_america", "west_asia")

 r <- raster(file.path(results, "africa", "BiomassFinalRaster.tif"))
 image(r, ylim = c(-90, 90), xlim = c(-180, 180))
 
for(reg in regions){

  r <- raster(file.path(results, reg, "BiomassFinalRaster.tif"))
  plot(r, add = TRUE)  
}

##############################################################################################
all_data <-"3_Data"
files <- list.files(file.path(all_data))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
all_data <- read.csv(file.path(all_data, paste("sites_",date,".csv", sep = "")))


models <- "Models"
load(file.path(models, "richnessmodel_full.rds"))
load(file.path(models, "biomassmodel_full.rds"))
load(file.path(models, "abundancemodel_full.rds"))

studies1 <- as.vector(unique(richness_model@frame$Study_Name))
studies2 <- as.vector(unique(abundance_model@frame$Study_Name))
studies3 <- as.vector(unique(biomass_model@frame$Study_Name))
all_studies <- c(studies1, studies2, studies3)
all_studies <- unique(all_studies)



all_studies <- all_data[all_data$Study_Name %in% all_studies,]

coord<-aggregate(cbind(all_studies$Longitude__Decimal_Degrees, all_studies$Latitude__decimal_degrees), list(all_studies$Study_Name), mean)

coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")



png(file.path(figures, "SpeciesRichness_DataPoints.png"),width=17.5,height=8.75,units="cm",res=1200)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))
image(spR_finalraster, col=b.cols(length(blubrks)-1), breaks=blubrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(spR_finalraster, col=r.cols(length(redbrks)-1), add = TRUE, breaks=redbrks, xaxt="n", yaxt="n", ylab="", xlab="")
points(dsSPDF, col="black", bg="black", cex= 0.6, pch=19)
#hist(1:10, axes = FALSE)
par(mar=c(1,13,1,13))
blu_scale <- rep(b.cols(length(blubrks)-1), each = 2) # 200 total
red_scale <-rep(r.cols(length(redbrks)-1), each = 2)



barplot(rep(1, 396), col = c(blu_scale, red_scale), border =c(blu_scale, red_scale), axes = FALSE )
mtext("1", at = 20, cex = 1)
# mtext("6", at = 350, cex = 0.5)
mtext("8", at = 455, cex = 1)

dev.off()


