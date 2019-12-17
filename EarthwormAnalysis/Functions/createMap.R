createSizedMap <- function(dat)
{
  
  try(if(!("Longitude__Decimal_Degrees" %in% names(dat))) stop("No longitude column - must be a specific name"))
  try(if(!("Latitude__decimal_degrees" %in% names(dat))) stop("No latitude column - must be a specific name"))
  
  if(any(is.na(dat$Latitude__decimal_degrees))){print("Some NAs in Latitude column")}
  if(any(is.na(dat$Longitude__Decimal_Degrees))){print("Some NAs in longitude column")}
  
  
  coord<-aggregate(cbind(dat$Longitude__Decimal_Degrees, dat$Latitude__decimal_degrees), list(dat$Study_Name), mean)
  coord$X<-coord$Group.1
  coord<-coord[2:4]
  names(coord)<-c("Long", "Lat", "X")
  
  dsSPDF<-SpatialPointsDataFrame(coord[,1:2], data.frame(coord[,1:3]))
  #proj4string(dsSPDF)<-CRS("+proj=longlat")
  proj4string(dsSPDF)<-CRS("+init=ESRI:54030")
  all(coord$X %in% names(table(dat$Study_Name)))
  studyN <- data.frame(table(dat$Study_Name))
  coord <- merge(coord, studyN, all.x = TRUE, by.x = "X", by.y = "Var1")
  
  names(coord)<-c("X", "Long", "Lat", "nSites")
  
  coord$size <- ((coord$nSites-min(coord$nSites))/(max(coord$nSites)-min(coord$nSites)) * 2) + 0.5
  
  dsSPDF<-SpatialPointsDataFrame(coord[,2:3], data.frame(coord[,1:5]))
  proj4string(dsSPDF)<-CRS("+proj=longlat")
  
  transpBlack <- rgb(0, 0, 0, alpha = 0.4, names = NULL, maxColorValue = 1)

  mar=c(0,0,0,0)
  map("world",border="gray87",fill=TRUE, col="gray87",mar=rep(0,4))
  points(dsSPDF, col=transpBlack, bg = transpBlack, cex= coord$size, pch=19)
  # corner.label2(label = "A", x = -1, y = 1, cex = plotlabcex, font = 2)
  
  sizes <- c(1, 50, 100, 150, 200, 250)
  cexsizes <- ((sizes-min(coord$nSites))/(max(coord$nSites)-min(coord$nSites)) * 2) + 0.5
  
  legend(x = -170, y = 2, legend = sizes, pch = 19, pt.cex =cexsizes, bty="n", cex = 0.7, 
         y.intersp = c(1, 1, 1, 1.05, 1.1, 1.18),
         x.intersp = c(1.19),
         title = "Number Of Sites")
  
}