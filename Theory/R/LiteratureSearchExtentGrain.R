setwd("C:/Users/hp39wasi/sWorm/Theory")

library(RColorBrewer)
library(scales)

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figs <- "Figures"



dat <- read.csv("LitSearch_data.csv")
table(dat$Grain.size..smallest.sample.unit.)
table(dat$Grain.size.unit)


table(dat$Extent)

### Fixing theory part
dat$Theory <- as.character(dat$Theory)
dat$Theory [which(dat$Theory == "IBT_SAR")] <- "IBT"
dat$Theory [which(dat$Theory == "MC")] <- "Meta"
dat$Theory [which(dat$Theory == "Neutral theory")] <- "Niche/Neutral"
dat$Theory [which(dat$Theory == "PROD")] <- "SER"

# "SER", "IBT", "MT", "SAR", "Niche/Neutral", "Meta"


dat <- droplevels(dat[!(dat$Theory %in% c("MT", "SAR")),])

## TODO
## There are four rows without extent, which I have now lost with the next line of code
dat$Extent[dat$Extent == ""] <- "unknown"
dat$Extent <- droplevels(dat$Extent)
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

dat$grain <- factor(dat$grain, levels = c("1cm", "10cm", "1m", "10m", "100m", "Unknown"))


dat$grainN <- as.numeric(dat$grain)
tapply(dat$grainN, dat$grain, mean)

dat$ExtentN <- as.numeric(dat$Extent)
tapply(dat$ExtentN, dat$Extent, mean)

plot(jitter(dat$grainN), jitter(dat$ExtentN))
## That works

### Now to add in sizes
levels(dat$Species.Group)
dat$Species <- dat$Species.Group

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

dat$Species <- droplevels(dat$Species)

dat$SizeGroup <- dat$Species
levels(dat$SizeGroup)[levels(dat$SizeGroup) %in% c("Earthworms", "soil invertebrates")] <- "Macrofauna"
levels(dat$SizeGroup)[levels(dat$SizeGroup) %in% c("mites", "springtails")] <- "Mesofauna"
levels(dat$SizeGroup)[levels(dat$SizeGroup) %in% c("nematodes", "Protozoa")] <- "Microfauna"
levels(dat$SizeGroup)[levels(dat$SizeGroup) %in% c("Fungi", "prokaryotes", "Soil microbes")] <- "Microorganisms"


# save file
write.csv(dat, file = "Data\\LitSearchFormatted.csv", row.names = FALSE)





#### With size groups coloured
#### Grain on y, extent on x
#### Without unknowns

dat <- droplevels(dat[dat$Extent != "unknown",])
dat <- droplevels(dat[dat$grain != "Unknown",])

colors <- brewer.pal(4, "Dark2")

dat$SizeGroup <- factor(dat$SizeGroup, levels= c("Microorganisms", "Microfauna", "Mesofauna", "Macrofauna"))

dat$SizeColours <- colors[1]

dat$SizeColours[dat$SizeGroup == "Microfauna"] <- colors[2]
dat$SizeColours[dat$SizeGroup == "Mesofauna"] <- colors[3]
dat$SizeColours[dat$SizeGroup == "Macrofauna"] <- colors[4]

pdf(file = file.path(figs, "GrainExtent_sizeClasses.pdf"), width = 11)
par(mar=c(4.5, 4, 1, 1))
layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
plot(jitter(grainN) ~ jitter(ExtentN), cex = 1.5,  pch = 20, col = SizeColours, data = dat,
     xaxt='n',  yaxt='n', xlab = "Extent", ylab = "Grain", bty = "l")
axis(2, at = 1:length(levels(dat$grain)), labels = levels(dat$grain), las = 2)
axis(1, at = 1:length(levels(dat$Extent)), labels = levels(dat$Extent), las = 2)
# Legend
par(mar=c(0,0,0,0))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 5), ylim = c(0, 1))
points(1:4, rep(0.5, 4) , cex = 1.5, col = colors, pch=19)
text(1:4, rep(0.5,4), labels = levels(dat$SizeGroup), pos =4)
dev.off()


# plot(jitter(grainN) ~ jitter(ExtentN), pch = 19, col = SizeColours, data = dat)

cols <-  brewer.pal(4, "Set1")
dat$TheoryColor <-  cols[1]
dat$TheoryColor[dat$Theory == "Meta"] <- cols[2]
dat$TheoryColor[dat$Theory == "Niche/Neutral"] <- cols[3]
dat$TheoryColor[dat$Theory == "SER"] <- cols[4]

dat$Theory <- as.factor(dat$Theory)

pdf(file = file.path(figs, "GrainExtent_theory.pdf"), width = 11)
par(mar=c(4.5, 4, 1, 1))
layout(matrix(c(1,2), 2,1, byrow = TRUE), c(5, 1), c(5, 1))
plot(jitter(grainN) ~ jitter(ExtentN), cex = 1.5,  pch = 20, col = TheoryColor, data = dat,
     xaxt='n',  yaxt='n', xlab = "Extent", ylab = "Grain", bty = "l")
axis(2, at = 1:length(levels(dat$grain)), labels = levels(dat$grain), las = 2)
axis(1, at = 1:length(levels(dat$Extent)), labels = levels(dat$Extent), las = 2)
# Legend
par(mar=c(0,0,0,0))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 5), ylim = c(0, 1))
points(1:4, rep(0.5, 4) , cex = 1.5, col = cols, pch=19)
text(1:4, rep(0.5,4), labels = levels(dat$Theory), pos =4)
dev.off()

########################
## When papers were published

all <- read.csv("AllSoilSearchResults.csv")

byyear <- as.data.frame(table(all$PY, all$Theory))
names(byyear) <- c("year", "theory", "freq")

library(ggplot2)
pdf(file = file.path(figs, "PublicationYear_theory.pdf"), width = 11)
ggplot(data = byyear, aes(x = factor(year), y = freq, color = theory)) +       
  geom_line(aes(group = theory)) + geom_point() +
  labs(x = "Year", y = "Frequency", colour = "")
dev.off()
#################
dat$Species <- factor(dat$Species, levels = c("prokaryotes","Fungi", "Soil microbes","nematodes",         
                                              "Protozoa", "mites","springtails", "soil invertebrates",
                                               "Earthworms"))

colors <- brewer.pal(11, "Spectral")
key <- data.frame(species = unique(dat$Species)[order(unique(dat$Species))], 
                  size = seq(0.5, 2.5, length = 9),
                  colour = colors[c(1:4, 7:11)])

# dat$speciesN <- as.numeric((dat$Species))
# tapply(dat$speciesN, dat$Species, mean)
# dat$speciesN <- (dat$speciesN + 10)/10

dat <- (merge(dat, key, by.x = "Species", by.y = "species"))

pdf(file = file.path(figs, "GrainExtent_alldata.pdf"), width = 11)
par(mar=c(4.5, 4, 1, 1), mfrow = c(1, 2))
plot(jitter(dat$grainN), jitter(dat$ExtentN), cex = dat$size,  #col = as.character(dat$colour), 
     xaxt='n',  yaxt='n', xlab = "Grain", ylab = "Extent", bty = "l")
axis(1, at = 1:length(levels(dat$grain)), labels = levels(dat$grain), las = 2)
axis(2, at = 1:length(levels(dat$Extent)), labels = levels(dat$Extent), las = 1)
# Legend
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 10), ylim = c(0, 10))
points(rep(1, 9), 1:9, cex = key$size)
text(rep(1.5, 9), 1:9, labels = key$species, pos =4)
dev.off()


dat_nounknown <- droplevels(dat[dat$grain != "Unknown",])
dat_nounknown <- droplevels(dat_nounknown[dat_nounknown$Extent != "Unknown",])

pdf(file = file.path(figs, "GrainExtent_no-unknowns.pdf"), width = 11)
par(mar=c(4.5, 4, 1, 1), mfrow = c(1, 2))
plot(jitter(dat_nounknown$grainN), jitter(dat_nounknown$ExtentN), cex = dat_nounknown$size,  #col = as.character(dat$colour), 
     xaxt='n',  yaxt='n', xlab = "Grain", ylab = "Extent", bty = "l")
axis(1, at = 1:length(levels(dat_nounknown$grain)), labels = levels(dat_nounknown$grain), las = 2)
axis(2, at = 1:length(levels(dat_nounknown$Extent)), labels = levels(dat_nounknown$Extent), las = 1)
# Legend
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 10), ylim = c(0, 10))
points(rep(1, 9), 1:9, cex = key$size)
text(rep(1.5, 9), 1:9, labels = key$species, pos =4)
dev.off()

###################
## With size on the Y
#####################

key_size <- data.frame(species = c("prokaryotes","Fungi","Algae","nematodes",
                                   "Protozoa","springtails","mites","Diplopoda","Earthworms"),
                         size_mm = c(0.001396368361,0.01445439771,0.01506607066,0.01936421964,
                                  0.02818382931,0.4246195639,0.5457578611,6.683439176,8.222426499))

dat_size <- droplevels(dat[!(dat$Species %in% c("Soil microbes",  "soil invertebrates")),])

dat_size <- merge(dat_size, key_size, by.x = "Species", by.y = "species")



key_theory <- data.frame(theory = c("SER", "IBT", "MT", "SAR", "Niche/Neutral", "Meta"),
                             colourTheory = brewer.pal(6, "Dark2")) 


dat_size <- merge(dat_size, key_theory, by.x = "Theory", by.y = "theory")

pdf(file = file.path(figs, "GrainExtent_bodysize.pdf"), width = 11)

par(mar=c(4.5, 4.5, 1, 1), mfrow = c(1, 3))
plot(jitter(dat_size$grainN), jitter(log10(dat_size$size_mm)), cex = 1.5,  
     col = alpha(as.character(dat_size$colourTheory),0.3) , pch =19,
     xaxt='n',  yaxt='n', xlab = "Grain", ylab = "log(Body Size, mm)", cex.lab = 1.5,  bty = "l")
axis(1, at = 1:length(levels(dat_size$grain)), labels = levels(dat_size$grain), las = 2)
axis(2, at = -3:1, labels = c("0.001", "0.01", "0.1", "1", "10"))

plot(jitter(dat_size$ExtentN), jitter(log10(dat_size$size_mm)), cex = 1.5,  
     col = alpha(as.character(dat_size$colourTheory),0.3) , pch =19,
     xaxt='n',  yaxt='n', xlab = "Extent", ylab = "log(Body Size, mm)", cex.lab = 1.5,  bty = "l")

axis(1, at = 1:length(levels(dat_size$Extent)), labels = levels(dat_size$Extent), las = 2)
axis(2, at = -3:1, labels = c("0.001", "0.01", "0.1", "1", "10"))

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 10), ylim = c(0, 10))
points(rep(1, 6), 1:6, cex = 1.5, pch = 19, col = alpha(as.character(key_theory$colourTheory),0.3))
text(rep(1.5, 9), 1:6, labels = key_theory$theory, pos =4)
dev.off()

