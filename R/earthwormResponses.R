## Directory with raw data
dir <- "C:/Users/hp39wasi/Dropbox/sWorm/survey/rawData"
outDir <- "C:/Users/hp39wasi/Dropbox/sWorm/survey/earthwormsOnly"

## Extracting most recent CSV
files <- list.files(dir)
file_dates <- sapply(strsplit(files, "_"), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates)
all <- read.csv(file.path(dir, paste(date, "_Soil biodiversity survey short-report.csv", sep="")))

## Removing all but those willing to share earthworm data
col <- "We.are.planning.a.detailed.analysis.of.global.earthworm.diversity.as.part.of.synthesis.group.at.sDiv.in.Germany..Do.you.have.earthworm.data.that.you.would.be.willing.to.contribute.to.this.project."
colnum <- which(names(all) == col)
all <- all[!(is.na(all[,colnum])),]
worms <- droplevels(all[all[,colnum] == 1,])

## Saving the results
write.csv(worms, file = file.path(outDir, paste(date, "_Soil biodiversity survey short-reportEarthworms.csv", sep="")))
