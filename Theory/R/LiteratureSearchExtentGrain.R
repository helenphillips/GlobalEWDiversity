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
dat$Theory [which(dat$Theory == "MC")] <- "MCT"
dat$Theory [which(dat$Theory == "Neutral theory")] <- "Niche/Neutral"
dat$Theory [which(dat$Theory == "PROD")] <- "SER"

# "SER", "IBT", "MT", "SAR", "Niche/Neutral", "Meta"


dat <- droplevels(dat[!(dat$Theory %in% c("MT", "SAR")),])
dat$Theory <- as.factor(dat$Theory)
dat$Theory <- factor(dat$Theory, levels = c("SER", "IBT", "MCT", "Niche/Neutral"))



dat$Extent <- as.character(dat$Extent)
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

dat$grain <- factor(dat$grain, levels = c("1cm", "10cm", "1m", "10m", "100m", "Unknown"))


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


pdf(file = file.path(figs, "MegaPlot.pdf"), width = 11)
m <- matrix(c(1, 2, 3, 4, 5, 6, 6, 6, 6, 7, 6, 6, 6, 6, 7), byrow = TRUE, ncol = 5)
layout(m)

par(mar=c(0.1, 0.1, 0.1, 0.1))
pie(sp$Freq[sp$Var2 == "SER"], labels = "", clockwise = TRUE, col=colors)
mtext("n = 15", 1)
pie(sp$Freq[sp$Var2 == "IBT"], labels = "", clockwise = TRUE, col=colors)
mtext("n = 15", 1)
pie(sp$Freq[sp$Var2 == "MCT"], labels = "", clockwise = TRUE, col=colors)
mtext("n = 18", 1)
pie(sp$Freq[sp$Var2 == "Niche/Neutral"], labels = "", clockwise = TRUE, col=colors)
mtext("n = 8", 1)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(1, 4.5), ylim = c(0,4))
points(rep(1, length = 4), 4:1, col = colors, pch = 19, cex = 2)
text(rep(1.2, length = 4), 4:1, adj = c(0, 0.5), labels = levels(dat$SizeGroup))


symb <- 1

par(mar=c(3, 4.5, 1, 4))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(1, 4.5), ylim = c(0, 6))
points(rep(1, length = nrow(x[x$Var2 == "SER",])), x$Var1[x$Var2 == "SER"], cex = log(x$Freq[x$Var2 == "SER"] + 1), pch = symb)
points(rep(2, length = nrow(x[x$Var2 == "IBT",])), x$Var1[x$Var2 == "IBT"], cex = log(x$Freq[x$Var2 == "IBT"] + 1), pch = symb)
points(rep(3, length = nrow(x[x$Var2 == "MCT",])), x$Var1[x$Var2 == "MCT"], cex = log(x$Freq[x$Var2 == "MCT"] + 1), pch = symb)
points(rep(4, length = nrow(x[x$Var2 == "Niche/Neutral",])), x$Var1[x$Var2 == "Niche/Neutral"], cex = log(x$Freq[x$Var2 == "Niche/Neutral"] + 1), pch = symb)

points(rep(1.4, length = nrow(y[y$Var2 == "SER",])), y$Var1[y$Var2 == "SER"], cex = log(y$Freq[y$Var2 == "SER"] + 1), pch = 19, col = "gray")
points(rep(2.4, length = nrow(y[y$Var2 == "IBT",])), y$Var1[y$Var2 == "IBT"], cex = log(y$Freq[y$Var2 == "IBT"] + 1), pch = 19, col = "gray")
points(rep(3.4, length = nrow(y[y$Var2 == "MCT",])), y$Var1[y$Var2 == "MCT"], cex = log(y$Freq[y$Var2 == "MCT"] + 1), pch = 19, col = "gray")
points(rep(4.4, length = nrow(y[y$Var2 == "Niche/Neutral",])), y$Var1[y$Var2 == "Niche/Neutral"], cex = log(y$Freq[y$Var2 == "Niche/Neutral"] + 1), pch = 19, col = "gray")

abline(v=c(1.7, 2.7, 3.7), lty = 2)

axis(1, at = c(1.2, 2.2, 3.2, 4.2), labels = c("SER", "IBT", "MCT", "Niche/Neutral"))
axis(2, at = 1:5, labels = levels(dat$grain), las = 2)
axis(4, at = 1:6, labels = levels(dat$Extent), las = 2)

mtext("Grain (Black circles)", side = 2, line = 3)
mtext("Extent (Gray dots)", side = 4, line = 3.5, col = "gray")

par(mar=c(0.1, 0.1, 0.1,0.1))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(1, 4.5), ylim = c(0,4))
points(rep(1.5, length = 4), c(3.5, 3, 2.5, 2), col = "black", pch = symb, cex = c(log(2), log(6), log(11), log(14)))
text(rep(1.8, length = 4), c(3.5, 3, 2.5, 2), adj = c(0, 0.5), labels = c(1, 5, 10, 13))
text(1.3, 3.8, "Number of Cases", adj = c(0, 0.5))

dev.off()



#############################################
## Grain/Extent on Y, Theory on X
#############################################

colors=c("#4477AA", "#117733", "#DDCC77", "#CC6677")

all$SizeGroup_color <- colors[1] ## microorganisms
all$SizeGroup_color[levels(all$SizeGroup) == "Microfauna"] <- colors[2]
all$SizeGroup_color[levels(all$SizeGroup) == "Mesofauna"] <- colors[3]
all$SizeGroup_color[levels(all$SizeGroup) == "Macrofauna"] <- colors[4]



ex <- data.frame(table(all$ExtentN, all$SizeGroup_color, all$Theory))
ex <- ex[ex$Freq > 0,]

gr <- data.frame(table(all$grainN, all$SizeGroup_color, all$Theory))
gr <- gr[gr$Freq > 0,]
gr <- gr[gr$Var1 != 6,] ## the unknowns

unknowns <- all[all$grainN == 6,]
u <- data.frame(table(unknowns$SizeGroup))

pdf(file = file.path(figs, "Extent+Grain+SizeClass.pdf"), width = 11)

layout(matrix(c(1, 0, 2, 3, 2, 3, 2, 3, 4,4), byrow = TRUE, ncol = 2))
par(mar=c(0, 0, 0, 0))
pie(u$Freq, labels = "", clockwise = TRUE, col=colors)
mtext("Number of unknown grain size = 20", 1)
par(mar=c(2, 4.5, 1, 1))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(1, 4), ylim = c(0,5))

points(jitter(rep(1, length = nrow(gr[gr$Var3 == "SER",]))), (gr$Var1[gr$Var3 == "SER"]), pch = 19, col = as.character(gr$Var2[gr$Var3 == "SER"]), cex = log(gr$Freq[gr$Var3 == "SER"] + 1))
points(jitter(rep(2, length = nrow(gr[gr$Var3 == "IBT",]))), (gr$Var1[gr$Var3 == "IBT"]), pch = 19, col = as.character(gr$Var2[gr$Var3 == "IBT"]), cex = log(gr$Freq[gr$Var3 == "IBT"] + 1))
points(jitter(rep(3, length = nrow(gr[gr$Var3 == "MCT",]))), (gr$Var1[gr$Var3 == "MCT"]), pch = 19, col = as.character(gr$Var2[gr$Var3 == "MCT"]), cex = log(gr$Freq[gr$Var3 == "MCT"] + 1))
points(jitter(rep(4, length = nrow(gr[gr$Var3 == "Niche/Neutral",]))), (gr$Var1[gr$Var3 == "Niche/Neutral"]), pch = 19, col = as.character(gr$Var2[gr$Var3 == "Niche/Neutral"]), cex = log(gr$Freq[gr$Var3 == "Niche/Neutral"] + 1))
axis(1, at = c(1, 2, 3, 4), labels = c("SER", "IBT", "MCT", "Niche/Neutral"))
axis(2, at = 1:5, labels = levels(dat$grain)[1:5], las = 2)
mtext("Grain", side = 2, line = 3.2)

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(1, 4), ylim = c(0,8))
points(jitter(rep(1, length = nrow(ex[ex$Var3 == "SER",]))), (ex$Var1[ex$Var3 == "SER"]), pch = 19, col = as.character(ex$Var2[ex$Var3 == "SER"]), cex = log(ex$Freq[ex$Var3 == "SER"] + 1))
points(jitter(rep(2, length = nrow(ex[ex$Var3 == "IBT",]))), (ex$Var1[ex$Var3 == "IBT"]), pch = 19, col = as.character(ex$Var2[ex$Var3 == "IBT"]), cex = log(ex$Freq[ex$Var3 == "IBT"] + 1))
points(jitter(rep(3, length = nrow(ex[ex$Var3 == "MCT",]))), (ex$Var1[ex$Var3 == "MCT"]), pch = 19, col = as.character(ex$Var2[ex$Var3 == "MCT"]), cex = log(ex$Freq[ex$Var3 == "MCT"] + 1))
points(jitter(rep(4, length = nrow(ex[ex$Var3 == "Niche/Neutral",]))), (ex$Var1[ex$Var3 == "Niche/Neutral"]), pch = 19, col = as.character(ex$Var2[ex$Var3 == "Niche/Neutral"]), cex = log(ex$Freq[ex$Var3 == "Niche/Neutral"] + 1))
axis(1, at = c(1, 2, 3, 4), labels = c("SER", "IBT", "MCT", "Niche/Neutral"))
axis(2, at = 1:8, labels = levels(all$Extent)[1:8], las = 2)
mtext("Extent", side = 2, line = 3.2)

# Legend
par(mar=c(0,0,0,0))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 5), ylim = c(0, 1))
points(1:4, rep(0.5, 4) , cex = 1.5, col = colors, pch=19)
text(1:4, rep(0.5,4), labels = levels(dat$SizeGroup), pos =4)
text(0.7, 0.5, labels= "Species Size Class")

points(1:4, rep(0.2, 4), col = "black", pch = 19, cex = c(log(2), log(3), log(4), log(5)))
text(1:4, rep(0.2, 4), pos = 4 , labels = c(1,2,3,4))
text(0.7, 0.23, "Number of Cases")


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
