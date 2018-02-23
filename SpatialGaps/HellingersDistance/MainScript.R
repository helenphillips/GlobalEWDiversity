################################
## DATA ##
##################################

## You would need to source the HellsDistFunction.R
## And load your data
## For this code, I have calcualted the bias between different land uses (primary, secondary, pasture and cropland)
## And for continuous variables of human population density, and distance to roads

## The code will also need the global layers
## In my functions "pri_1km_int" etc. are the global layers.




sites$IntDistToRoads <- round(sites$DistToRoad)

primary <- droplevels(sites[sites$Predominant_habitat == "Primary vegetation",])
secondary <- droplevels(sites[sites$Predominant_habitat == "Secondary vegetation",])
pasture <- droplevels(sites[sites$Predominant_habitat == "Pasture",])
cropland <- droplevels(sites[sites$Predominant_habitat == "Cropland",])
urban <- droplevels(sites[sites$Predominant_habitat == "Urban",])

################################
## LIBRARIES ##
##################################

library(ape)
library(distrEx)
library(rgdal)
library(ncf)
library(Hmisc)
library(gplots)
#################################
### ANALYSIS
##################################


################################################
## Method used in Gonzales et al 2016
quants <- c(1:4)
UD<- DiscreteDistribution(quants, c(0.25,0.25,0.25,0.25))

pri_HD <- list()
pri_HD[[1]] <- HellsDist(site_data_col = primary$Primary_prop, "pri_1km_int", inDir=inDir, land.use = TRUE, n.reps = 1000, UD = UD)
pri_HD[[2]] <- HellsDistData(site_data_col = primary$Primary_prop,  "pri_1km_int", inDir = inDir, land.use = TRUE, UD = UD)

sec_HD <- list()
sec_HD[[1]] <- HellsDist(site_data_col = secondary$Secondary_prop, "sec_1km_int", inDir=inDir, land.use = TRUE, n.reps = 1000, UD = UD)
sec_HD[[2]] <- HellsDistData(site_data_col = secondary$Secondary_prop,  "sec_1km_int", inDir = inDir, land.use = TRUE, UD = UD)

pas_HD <- list()
pas_HD[[1]] <- HellsDist(site_data_col = pasture$Pasture_prop, "pas_1km_int", inDir=inDir, land.use = TRUE, n.reps = 1000, UD = UD)
pas_HD[[2]] <- HellsDistData(site_data_col = pasture$Pasture_prop,  "pas_1km_int", inDir = inDir, land.use = TRUE, UD = UD)

crp_HD <- list()
crp_HD[[1]] <- HellsDist(site_data_col = cropland$Crop_prop, "crp_1km_int", inDir=inDir, land.use = TRUE, n.reps = 1000, UD = UD)
crp_HD[[2]] <- HellsDistData(site_data_col = cropland$Crop_prop,  "crp_1km_int", inDir = inDir, land.use = TRUE, UD = UD)

urb_HD <- list()
urb_HD[[1]] <- HellsDist(site_data_col = urban$Urban_prop, "urb_1km_int", inDir=inDir, land.use = TRUE, n.reps = 1000, UD = UD)
urb_HD[[2]] <- HellsDistData(site_data_col = urban$Urban_prop,  "urb_1km_int", inDir = inDir, land.use = TRUE, UD = UD)

global_HD <- list()
global_HD[[1]] <- HellsDist(site_data_col = sites$HPD, "hpd", inDir="C:\\Helens\\CountrysideBiogeography\\RasterOutputs\\HPD\\", land.use = FALSE, n.reps = 1000, UD = UD)
global_HD[[2]] <- HellsDistData(site_data_col = sites$HPD,  "hpd", inDir="C:\\Helens\\CountrysideBiogeography\\RasterOutputs\\HPD\\", land.use = FALSE, UD = UD)
global_HD[[3]] <- HellsDist(site_data_col = sites$IntDistToRoads, "Int-rddistwgs", inDir="C:\\Helens\\CountrysideBiogeography\\RasterOutputs\\DistToRoads\\", land.use = FALSE, n.reps = 1000, UD = UD)
global_HD[[4]] <- HellsDistData(site_data_col = sites$IntDistToRoads,  "Int-rddistwgs", inDir="C:\\Helens\\CountrysideBiogeography\\RasterOutputs\\DistToRoads\\", land.use = FALSE, UD = UD)

save(list=c("crp_HD", "pas_HD", "pri_HD", "sec_HD", "urb_HD", "global_HD"), file=paste(path, "Results\\HellingersDist.rda", sep = ""))
load(file.path(path, "Results", "HellingersDist.rda"))


#####################
###### Plots
#####################

pri_HD_meansSD <- c(mean(pri_HD[[1]]), sd(pri_HD[[1]]))
pridat_HD_meansSD <- c(pri_HD[[2]])

sec_HD_meansSD <- c(mean(sec_HD[[1]]), sd(sec_HD[[1]]))
secdat_HD_meansSD <- c(sec_HD[[2]])

pas_HD_meansSD <- c(mean(pas_HD[[1]]), sd(pas_HD[[1]]))
pasdat_HD_meansSD <- c(pas_HD[[2]])

crp_HD_meansSD <- c(mean(crp_HD[[1]]), sd(crp_HD[[1]]))
crpdat_HD_meansSD <- c(crp_HD[[2]])

urb_HD_meansSD <- c(mean(urb_HD[[1]]), sd(urb_HD[[1]]))
urbdat_HD_meansSD <- c(urb_HD[[2]])

global_HD_meansSD <- c(mean(global_HD[[1]]), mean(global_HD[[3]]), sd(global_HD[[1]]), sd(global_HD[[3]]))
globaldat_HD_meansSD <- c(global_HD[[2]], global_HD[[4]])

pdf(file = file.path(path, "Figures", "HellDist.pdf"))
par(mar=c(3, 7, 1, 1))
pt.cex <- 0.5
pt.pch <- 19
layer_labels <- c("Primary", "Secondary", "Pasture", "Cropland", "Urban", "HPD", "DistToRoads")
plotCI(x=1, y=0, type="n", xlim=c(0,0.61), ylim=c(1,7), xlab = "Hellinger's d", yaxt = "n", ylab = "")
rect(-0.5,2.5, 1, 8,col="#A0A0A033",border=NA)

plotCI(x =  pri_HD_meansSD[1], y =7, uiw = pri_HD_meansSD[2],err='x', cex = pt.cex, gap =0, pch = pt.pch, add =TRUE)
plotCI(x =  pridat_HD_meansSD[1], y =7,err='x',cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)

plotCI(x =  sec_HD_meansSD[1], y =6, uiw = sec_HD_meansSD[2],err='x', cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)
plotCI(x =  secdat_HD_meansSD[1], y =6,err='x',cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)

plotCI(x =  pas_HD_meansSD[1], y =5, uiw = pas_HD_meansSD[2],err='x', cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)
plotCI(x =  pasdat_HD_meansSD[1], y =5,err='x',cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)

plotCI(x =  crp_HD_meansSD[1], y =4, uiw = crp_HD_meansSD[2],err='x', cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)
plotCI(x =  crpdat_HD_meansSD[1], y =4,err='x',cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)

plotCI(x =  urb_HD_meansSD[1], y =3, uiw = urb_HD_meansSD[2],err='x', cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)
plotCI(x =  urbdat_HD_meansSD[1], y =3,err='x',cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)

plotCI(x =  global_HD_meansSD[1:2], y =c(2, 1), uiw =  global_HD_meansSD[3:4],err='x', cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)
plotCI(x =  globaldat_HD_meansSD[1:2], y =c(2, 1),err='x',cex = pt.cex, gap =0, pch = pt.pch, add = TRUE)
axis(2, at = 7:1, labels = layer_labels, las =2)
dev.off()
#### At the end
(HD-(mean(HD.vals)))/(sd(HD.vals)) 