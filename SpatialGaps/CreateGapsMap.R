library(raster)
library(dismo)
library(maps)

sites <- read.csv("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\4_Data\\sitesRichness_2017-12-04.csv")

map <- raster("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\6_Data\\spRFinalRaster.tif")
mask <- raster("AccClimateMask.tif")

plot(map)



coord<-aggregate(cbind(sites$Longitude__Decimal_Degrees, sites$Latitude__decimal_degrees), list(sites$Study_Name), mean)
coord <- coord[complete.cases(coord),]
coord$X<-coord$Group.1
coord<-coord[2:4]
names(coord)<-c("Long", "Lat", "X")

dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")


cs <- circles(coord[, c(1,2)], d=100000, lonlat = TRUE)
plot(cs, add = TRUE)

p <- polygons(cs)
plot(p, add = TRUE, col = "black")



### Cities
library(rgdal)
cities <- readOGR(dsn = "C:\\Users\\hp39wasi\\Desktop\\World_Cities\\v10\\cities.gdb")

cities$rank <- as.factor(cities$POP_CLASS)
levels(cities$rank)[levels(cities$rank) == "Less than 50,000" ] <- "1"
levels(cities$rank)[levels(cities$rank) == "50,000 to 99,999"] <- "2"
levels(cities$rank)[levels(cities$rank) == "100,000 to 249,999"] <- "3"
levels(cities$rank)[levels(cities$rank) == "250,000 to 499,999" ] <- "4"
levels(cities$rank)[levels(cities$rank) == "500,000 to 999,999" ] <- "5"
levels(cities$rank)[levels(cities$rank) == "1,000,000 to 4,999,999" ] <- "6"
levels(cities$rank)[levels(cities$rank) == "5,000,000 and greater" ] <- "7"
cities$rank <- as.integer(as.character(cities$rank))


library("rgeos")
cities_buf <- gBuffer(cities, width = 100000, quadsegs = 10)

plot(cities, add = TRUE, pch = 19, cex = (cities$rank/10), col = "red")



###################################
## Actual Species Richness Map
###################################


library(RColorBrewer)

blubrks <- seq( -0.5,1, length.out = 100)
blues <- rev(brewer.pal(9, "Blues"))
b.cols <- colorRampPalette(blues, space = "rgb")


redbrks <-  seq(1.01, 2, length.out = 100)
reds <- brewer.pal(9, "Reds")
r.cols <- colorRampPalette(reds, space = "rgb")

png("SpeciesRichness_Gaps.png",width=17.5,height=8.75,units="cm",res=1200)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))
image(map, col=b.cols(length(blubrks)-1), breaks=blubrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(map, col=r.cols(length(redbrks)-1), add = TRUE, breaks=redbrks, xaxt="n", yaxt="n", ylab="", xlab="")



plot(cities, add = TRUE, pch = 19, cex = (cities$rank/10), col = "green")


plot(p, add = TRUE, col = "black")


#hist(1:10, axes = FALSE)
par(mar=c(1,13,1,13))
blu_scale <- rep(b.cols(length(blubrks)-1), each = 2) # 200 total
red_scale <-rep(r.cols(length(redbrks)-1), each = 2)



barplot(rep(1, 396), col = c(blu_scale, red_scale), border =c(blu_scale, red_scale), axes = FALSE )
mtext("1", at = 58, cex = 0.5)
# mtext("6", at = 350, cex = 0.5)
mtext("7", at = 408, cex = 0.5)

dev.off()

rm(sr)
gc()



############################################
## Change in climate versus sampled points

pdf(file = "ClimateChange_gaps.pdf")
map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
plot(mask, add = TRUE, legend = FALSE)
points(dsSPDF, col="black", bg="black", cex= 1, pch=19)
dev.off()
