setwd("C:/Users/hp39wasi/sWorm/Theory")


dat <- read.csv("Data\\LitSearchFormatted.csv")
dat <- dat[dat$Supported != "doesn't mention any, could just delete I think",]
dat <- droplevels(dat)
length(unique(dat$Title))

ds <- duplicated(dat$Title)
dat <- dat[!(ds),] ## non-duplicated titles

bib <- read.csv("AllSoilSearchResults.csv")




allDat <- merge(dat, bib, by.x = "Title", by.y = "TI", all.x = TRUE)
ds <- duplicated(allDat$Title)
allDat <- allDat[!(ds),]



names(allDat) ## I'll manually remove columns I think


write.csv(allDat, file = "Data\\Data+Bib.csv", row.names = FALSE)
