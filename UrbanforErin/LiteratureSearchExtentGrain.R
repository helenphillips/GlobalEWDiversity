# setwd("C:/Users/hp39wasi/sWorm/Theory")
setwd("C:/Users/hp39wasi/sWorm/UrbanforErin")
library(RColorBrewer)
library(scales)

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figs <- "Figures"



dat <- read.csv("../Theory/LitSearch_data.csv")
table(dat$Grain.size..smallest.sample.unit.)
table(dat$Grain.size.unit)


table(dat$Extent)

### Fixing theory part
dat$Theory <- as.character(dat$Theory)
dat$Theory [which(dat$Theory == "IBT_SAR")] <- "TIB"
dat$Theory [which(dat$Theory == "IBT")] <- "TIB"
dat$Theory [which(dat$Theory == "MC")] <- "MCT"
dat$Theory [which(dat$Theory == "Neutral theory")] <- "Niche/Neutral"
dat$Theory [which(dat$Theory == "PROD")] <- "SER"

dat <- droplevels(dat)
# "SER", "IBT", "MT", "SAR", "Niche/Neutral", "Meta"


dat <- droplevels(dat[!(dat$Theory %in% c("MT", "SAR")),])
dat$Theory <- as.factor(dat$Theory)
dat$Theory <- factor(dat$Theory, levels = c("SER", "TIB", "MCT", "Niche/Neutral"))



dat$Extent <- as.character(dat$Extent)
dat$Extent[dat$Extent == ""] <- "unknown"
dat$Extent[is.na(dat$Extent)] <- "unknown"
dat$Extent <- as.factor(dat$Extent)
# dat$Extent <- droplevels(dat$Extent)
dat$Extent <- factor(dat$Extent, levels = c("1m", "10m", "100m", "1km", "10km",
                                            "100km", "1000km", "Global", "unknown"))

dat$grain <- as.factor(paste(dat$Grain.size..smallest.sample.unit., dat$Grain.size.unit, sep =" "))
unique(dat$grain)
levels(dat$grain)[levels(dat$grain) == "unknown unknown"] <- "Unknown"  
levels(dat$grain)[levels(dat$grain) == "25 m2"] <- "10m"
levels(dat$grain)[levels(dat$grain) == "0.07 m"] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "2500 m2"] <- "100m"
levels(dat$grain)[levels(dat$grain) ==  "39 cm2"] <- "10cm"   
levels(dat$grain)[levels(dat$grain) == "1 litre" ]  <- "10cm"        
levels(dat$grain)[levels(dat$grain) == "0.15 ha" ] <- "10m" 
levels(dat$grain)[levels(dat$grain) == "0.01 m2"]  <- "10cm"    
levels(dat$grain)[levels(dat$grain) == "14 km2" ] <- "Unknown" ## This shouldn't be there anyway
levels(dat$grain)[levels(dat$grain) == "NA NA"] <- "Unknown"
levels(dat$grain)[levels(dat$grain) == "15 cm"] <- "10cm"
levels(dat$grain)[levels(dat$grain) =="3 cm"] <- "1cm"
levels(dat$grain)[levels(dat$grain) =="225 cm2"] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "seedling "] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "5 cm2"] <- "1cm"
levels(dat$grain)[levels(dat$grain) == "490 mm2"] <- "1cm"
levels(dat$grain)[levels(dat$grain) == "1 liter"] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "1 m2"] <- "10cm"
levels(dat$grain)[levels(dat$grain) =="25 cm2"] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "10 cm"] <- "10cm"
levels(dat$grain)[levels(dat$grain) =="0.25 m2"] <- "1m"           
levels(dat$grain)[levels(dat$grain) == "50 g"] <- "10cm"
levels(dat$grain)[levels(dat$grain) =="5 cm"] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "100 g"] <- "10cm"           
levels(dat$grain)[levels(dat$grain) =="40 mm diameter mm"] <- "1cm"
levels(dat$grain)[levels(dat$grain) == "2 cm diameter cm"] <- "1cm"
levels(dat$grain)[levels(dat$grain) == "1 g"] <- "1cm"
levels(dat$grain)[levels(dat$grain) == "10 g"] <- "1cm"
levels(dat$grain)[levels(dat$grain) == "100 mg"] <- "1cm"
levels(dat$grain)[levels(dat$grain) == "10 cm deep "] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "100 cm2"] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "individual leaves "] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "10 cm3"] <- "1cm"
levels(dat$grain)[levels(dat$grain) ==  "unknown "]  <- "Unknown"
levels(dat$grain)[levels(dat$grain) =="25 cm3"] <- "1cm"        
levels(dat$grain)[levels(dat$grain) == "20 cm2"] <-  "1cm"
levels(dat$grain)[levels(dat$grain) == "40 cm2"] <- "10cm"            
levels(dat$grain)[levels(dat$grain) == "16 cm2"] <- "1cm"
levels(dat$grain)[levels(dat$grain) == "625 cm2"] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "512 cm3"] <- "10cm"
levels(dat$grain)[levels(dat$grain) == "15625 cm3"] <- "10cm"    

dat$grain <- droplevels(dat$grain)

dat$grain <- factor(dat$grain, levels = c("1cm", "10cm", "1m", "10m", "Unknown"))


dat$grainN <- as.numeric(dat$grain)
tapply(dat$grainN, dat$grain, mean)

dat$ExtentN <- as.numeric(dat$Extent)
tapply(dat$ExtentN, dat$Extent, mean)

plot(jitter(dat$grainN), jitter(dat$ExtentN))
## That works

### Now to add in sizes
levels(dat$Species.Group)
dat$Species <- as.character(dat$Species.Group)

dat$Species[grep("fungi", dat$Species, ignore.case = TRUE, value = FALSE)] <- "Fungi"
dat$Species[grep("mycorrhiza", dat$Species, ignore.case = TRUE, value = FALSE)] <- "Fungi"
dat$Species[grep("Ascomycota", dat$Species, ignore.case = TRUE, value = FALSE)] <- "Fungi"
dat$Species[grep("bacteria", dat$Species, ignore.case = TRUE, value = FALSE)] <- "prokaryotes"
dat$Species[grep("Gram", dat$Species, ignore.case = TRUE, value = FALSE)] <- "prokaryotes"
dat$Species[grep("Archaea", dat$Species, ignore.case = TRUE, value = FALSE)] <- "prokaryotes"
dat$Species[grep("diazotrophic", dat$Species, ignore.case = TRUE, value = FALSE)] <- "prokaryotes"

dat$Species[grep("collembola", dat$Species, ignore.case = TRUE)] <- "springtails"
dat$Species[grep("springtails", dat$Species, ignore.case = TRUE, value = FALSE)] <- "springtails"
dat$Species[grep("earthworms", dat$Species, ignore.case = TRUE, value = FALSE)] <- "Earthworms"
dat$Species[grep("oribatid", dat$Species, ignore.case = TRUE, value = FALSE)] <- "mites"
dat$Species[grep("oribatida", dat$Species, ignore.case = TRUE, value = FALSE)] <- "mites"
dat$Species[grep("mesostigmatids", dat$Species, ignore.case = TRUE, value = FALSE)] <- "mites"
dat$Species[grep("amoebae", dat$Species, ignore.case = TRUE, value = FALSE)] <- "Protozoa"
dat$Species[grep("Microarthropods", dat$Species, ignore.case = TRUE, value = FALSE)] <- "soil invertebrates"
dat$Species[grep("Macroinvertebrates", dat$Species, ignore.case = TRUE, value = FALSE)] <- "soil invertebrates"
dat$Species[grep("litter arthropods", dat$Species, ignore.case = TRUE, value = FALSE)] <- "soil invertebrates"
dat$Species[grep("N Fixers", dat$Species, ignore.case = TRUE, value = FALSE)] <- "Soil microbes"



# dat$Species <- droplevels(dat$Species)

dat$SizeGroup <- as.factor(dat$Species)
levels(dat$SizeGroup)[levels(dat$SizeGroup) %in% c("Earthworms", "soil invertebrates")] <- "Macrofauna"
levels(dat$SizeGroup)[levels(dat$SizeGroup) %in% c("mites", "springtails")] <- "Mesofauna"
levels(dat$SizeGroup)[levels(dat$SizeGroup) %in% c("nematodes", "Protozoa")] <- "Microfauna"
levels(dat$SizeGroup)[levels(dat$SizeGroup) %in% c("Fungi", "prokaryotes", "Soil microbes")] <- "Microorganisms"

dat$SizeGroup <- factor(dat$SizeGroup, levels= c("Microorganisms", "Microfauna", "Mesofauna", "Macrofauna"))


# save file
write.csv(dat, file = "Data\\LitSearchFormatted.csv", row.names = FALSE)

all <- dat

#########################
## An alternative plot
#########################
dat <- droplevels(dat[dat$Extent != "unknown",])
dat <- droplevels(dat[dat$grain != "Unknown",])

# colors <- brewer.pal(11, "RdYlBu")
# colors <- colors[11:8]
colors=c("#4477AA", "#117733", "#DDCC77", "#CC6677")


tapply(dat$grainN, dat$Theory, mean)
tapply(dat$grainN, dat$Theory, summary)
x <- data.frame(table(dat$grainN, dat$Theory))
x <- x[x$Freq > 0,]

tapply(dat$ExtentN, dat$Theory, mean)
y <- data.frame(table(dat$ExtentN, dat$Theory))
y <- y[y$Freq > 0,]


sp <- data.frame(table(dat$SizeGroup, dat$Theory))

theories <- levels(dat$Theory)


# pdf(file = file.path(figs, "MegaPlot.pdf"), width = 11)
jpeg(file = file.path(figs, "MegaPlot.jpg"), width = 13.5, height = 10, units = "in", res = 300)

m <- matrix(c(1, 2, 3, 4, 5, 6, 6, 6, 6, 7, 6, 6, 6, 6, 7), byrow = TRUE, ncol = 5)
layout(m)

par(mar=c(0.1, 0.1, 0.1, 0.1))
pie(sp$Freq[sp$Var2 == theories[1]], labels = "", clockwise = TRUE, col=colors)
mtext(paste("n = ", nrow(dat[dat$Theory == theories[1],]), sep =""), 1, cex = 1.2)
pie(sp$Freq[sp$Var2 == theories[2]], labels = "", clockwise = TRUE, col=colors)
mtext(paste("n = ", nrow(dat[dat$Theory == theories[2],]), sep =""), 1, cex = 1.2)
pie(sp$Freq[sp$Var2 == theories[3]], labels = "", clockwise = TRUE, col=colors)
mtext(paste("n = ", nrow(dat[dat$Theory == theories[3],]), sep =""), 1, cex = 1.2)
pie(sp$Freq[sp$Var2 == theories[4]], labels = "", clockwise = TRUE, col=colors)
mtext(paste("n = ", nrow(dat[dat$Theory == theories[4],]), sep =""), 1, cex = 1.2)

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(1, 4.5), ylim = c(0,4))
points(rep(1, length = 4), 4:1, col = colors, pch = 19, cex = 3)
text(rep(1.2, length = 4), 4:1, adj = c(0, 0.5), labels = levels(dat$SizeGroup), cex = 2)


symb <- 1
par(mar=c(3, 7, 1, 8))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(1, 4.5), ylim = c(0, 6))
points(rep(1, length = nrow(x[x$Var2 == theories[1],])), x$Var1[x$Var2 == theories[1]], cex = log(x$Freq[x$Var2 == theories[1]] + 1.5), pch = symb)
points(rep(2, length = nrow(x[x$Var2 == theories[2],])), x$Var1[x$Var2 == theories[2]], cex = log(x$Freq[x$Var2 == theories[2]] + 1.5), pch = symb)
points(rep(3, length = nrow(x[x$Var2 == theories[3],])), x$Var1[x$Var2 == theories[3]], cex = log(x$Freq[x$Var2 == theories[3]] + 1.5), pch = symb)
points(rep(4, length = nrow(x[x$Var2 == theories[4],])), x$Var1[x$Var2 == theories[4]], cex = log(x$Freq[x$Var2 == theories[4]] + 1.5), pch = symb)

points(rep(1.4, length = nrow(y[y$Var2 == theories[1],])), y$Var1[y$Var2 == theories[1]], cex = log(y$Freq[y$Var2 == theories[1]] + 1.5), pch = 19)
points(rep(2.4, length = nrow(y[y$Var2 == theories[2],])), y$Var1[y$Var2 == theories[2]], cex = log(y$Freq[y$Var2 == theories[2]] + 1.5), pch = 19)
points(rep(3.4, length = nrow(y[y$Var2 == theories[3],])), y$Var1[y$Var2 == theories[3]], cex = log(y$Freq[y$Var2 == theories[3]] + 1.5), pch = 19)
points(rep(4.4, length = nrow(y[y$Var2 == theories[4],])), y$Var1[y$Var2 == theories[4]], cex = log(y$Freq[y$Var2 == theories[4]] + 1.5), pch = 19)

abline(v=c(1.7, 2.7, 3.7), lty = 2)

axis(1, at = c(1.2, 2.2, 3.2, 4.2), labels = theories, cex.axis =2)
axis(2, at = 1:length(unique(dat$grainN)), labels = levels(dat$grain), las = 2, cex.axis=2)
axis(4, at = 1:length(unique(dat$ExtentN)), labels = levels(dat$Extent), las = 2, cex.axis=2)

mtext("Grain (Black circles)", side = 2, line = 4.9, cex = 2)
mtext("Extent (Black dots)", side = 4, line = 7, col = "gray", cex = 2)

par(mar=c(0.1, 0.1, 0.1,0.1))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(1, 4.5), ylim = c(0,4))
points(rep(1.5, length = 4), c(3.5, 3, 2.5, 2), col = "black", pch = symb, cex = c(log(2.5), log(6.5), log(11.5), log(14.5)))
text(rep(1.8, length = 4), c(3.5, 3, 2.5, 2), adj = c(0, 0.5), labels = c(1, 5, 10, 13), cex =1.8)
text(1.3, 3.8, "Number of Cases", adj = c(0, 0.5), cex =2.2)

dev.off()

