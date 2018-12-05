setwd("C:/Users/hp39wasi/sWorm/Theory")

library(RColorBrewer)
library(scales)

if(!dir.exists("Figures")){
  dir.create("Figures")
}
figs <- "Figures"


dat <- read.csv("Data\\LitSearchFormatted.csv")

colors=c("#4477AA", "#117733", "#DDCC77", "#CC6677")


####################################################
# FUNCTIONS FOR FIGURES
####################################################

## PIECHART

pieChart <- function(dat, support = "Yes"){
  print("Function is expecting the full dataset for a theory")
  sp <- data.frame(table(dat$SizeGroup, dat$Support))
  
  pie(sp$Freq[sp$Var2 == support], labels = "", clockwise = TRUE, col=colors)
  mtext(paste("n = ", nrow(dat[dat$Supported == support,]), sep =""), 1)
}


ExtentAndGrain <- function(dat, support = "Yes"){
 
  grainLabs <- c("1cm", "10cm", "1m", "10m", "Unknown")
  extentLabs <- c("1m", "10m", "100m", "1km", "10km",
                  "100km", "1000km", "Global", "unknown")
  x <- data.frame(table(dat$grainN, dat$Supported))
  x <- x[x$Var2 == support,]
  x$Var2 <- NULL
  x <- x[x$Freq > 0,]
  new <- data.frame(Var1 = as.factor(setdiff(1:5, x$Var1)))
  new$Freq <- 0
  x <- rbind(x, new)
  x$Var1 <- as.integer(as.character(x$Var1))
  plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(1, 1.5), ylim = c(0, 9))
  points(rep(1, 5), x$Var1, cex = log(x$Freq + 1), pch = 1)
  axis(2, at = 1:length(grainLabs), labels = grainLabs, las = 2)
  
  y <- data.frame(table(dat$ExtentN, dat$Support))
  y <- y[y$Var2 == support,]
  y <- y[y$Freq > 0,]
  y$Var2 <- NULL
  new <- data.frame(Var1 = as.factor(setdiff(1:9, y$Var1)))
  new$Freq <- 0
  y <- rbind(y, new)
  y$Var1 <- as.integer(as.character(y$Var1))
  points(rep(1.4, length = nrow(y)), y$Var1, cex = log(y$Freq + 1), pch = 19)
  axis(4, at = 1:length(extentLabs), labels = extentLabs, las = 2)
  
}
########################################################
levels(dat$Supported)

dat$Supported[grep("yes - s", dat$Supported, ignore.case = TRUE)] <- "Yes"
dat$Supported[grep("yes", dat$Supported)] <- "Yes"
dat$Supported[grep("yes - m", dat$Supported, ignore.case = TRUE)] <- "Yes"

dat <- dat[dat$Supported != "doesn't mention any, could just delete I think",]

dat$Supported <- droplevels(dat$Supported)


jpeg(file.path(figs, "SupportOrNot.jpeg"), width = 3000, height = 2000, pointsize = 30)
## IBT
n <- layout(matrix(c(1,2,2,3,3,seq(4, 23, by = 1)), byrow = TRUE, ncol = 5))

par(mar=c(0,0,0,0))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, labels = "Theory")
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, labels = "Supported")
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, labels = "Not Support")

ibt <- droplevels(dat[dat$Theory == "TIB",])
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, labels = "IBT")
pieChart(ibt, support = "Yes")
ExtentAndGrain(ibt, support = "Yes")
pieChart(ibt, support = "No")
ExtentAndGrain(ibt, support = "No")

## SER
ser <- droplevels(dat[dat$Theory == "SER",])
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, labels = "SER")
pieChart(ser, support = "Yes")
ExtentAndGrain(ser, support = "Yes")
pieChart(ser, support = "No")
ExtentAndGrain(ser, support = "No")

## MCT
mct <- droplevels(dat[dat$Theory == "MCT",])
levels(mct$Supported)[levels(mct$Supported) == "No (all tested)"] <- "No"
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, labels = "MCT")

pieChart(mct, support = "Yes")
ExtentAndGrain(mct, support = "Yes")
pieChart(mct, support = "No")
ExtentAndGrain(mct, support = "No")

dev.off()


### If a yes for neutral, then it's a no for niche
### If a yes for niche, then its a no for neutral
# if its a no for either, we do not assume support for the other
