## Create map

spR_finalraster <- raster(file.path(soil_GLs, "spRFinalRaster.tif"))


library(RColorBrewer)

summary(spR_finalraster)


blubrks <- seq( -0.5,1, length.out = 100)
blues <- rev(brewer.pal(9, "Blues"))
b.cols <- colorRampPalette(blues, space = "rgb")


redbrks <-  seq(1.01, 2, length.out = 100)
reds <- brewer.pal(9, "Reds")
r.cols <- colorRampPalette(reds, space = "rgb")

png("SpeciesRichness.png",width=17.5,height=8.75,units="cm",res=1200)
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
mtext("1", at = 58, cex = 0.5)
# mtext("6", at = 350, cex = 0.5)
mtext("7", at = 408, cex = 0.5)

dev.off()

rm(sr)
gc()


####### Total Abundance Plot #######
ta <-readGDAL(paste(dir,"tabundance",sep=""))


pdf("C://Helens//CountrysideBiogeography//Results//Figures//TAfrequency.pdf")
par(mfrow=c(1, 2))
hist(ta$band1, xlab = "% Change", main="")
mtext("(a)", at = -50)
hist(ta$band1, xlab = "% Change", main="", ylim = c(0, 1500))
mtext("(b)", at = -50)
dev.off()


## Determining where to break the data
hist(ta$band1, xlab = "% Change", main="", ylim = c(0, 1e+07), xlim=c(-100, 0))
hist(ta$band1, xlab = "% Change", main="", ylim = c(0, 1e+03), xlim=c(0, 600))



blubrks <- c(seq( 1,100, length.out = 100), max(ta$band1, na.rm = TRUE))
blues <- brewer.pal(9, "Blues")
b.cols <- colorRampPalette(blues, space = "rgb")


redbrks <- c(min(ta$band1, na.rm = TRUE), seq(-20, -1, length.out = 100)) # This now groups everything below -20 into one bin. I think. 
reds <- rev(brewer.pal(9, "Reds"))
r.cols <- colorRampPalette(reds, space = "rgb")

png(paste(outDir,"TotalAbundance.png",sep=""),width=17.5,height=8.75,units="cm",res=1200)
nf <- layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
layout.show(nf)
par(mar=c(0.1,0.1,0.1,0.1))
image(ta, col=b.cols(length(blubrks)-1), breaks=blubrks, xaxt="n", yaxt="n", ylab="", xlab="")
image(ta, col=r.cols(length(redbrks)-1), add = TRUE, breaks=redbrks, xaxt="n", yaxt="n", ylab="", xlab="")
#hist(1:10, axes = FALSE)
par(mar=c(1,13,1,13))
blu_scale <- c(b.cols(length(blubrks)-1), rep(b.cols(length(blubrks)-1)[100], times = 75)) # 200 total
red_scale <-c(rep(r.cols(length(redbrks)-1)[1], times = 75), r.cols(length(redbrks)-1))

white <- rep("#FFFFFF", times = 2)

barplot(rep(1, 352), col = c(red_scale, white, blu_scale), border =c(red_scale, white, blu_scale), axes = FALSE )
mtext("-20%", at = 85, cex = 0.5)
mtext("100%", at = 330, cex = 0.5)
mtext("0%", at = 215, cex = 0.5)

dev.off()