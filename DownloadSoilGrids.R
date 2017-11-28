library(XML)
library(RCurl)

savedir <- "C:\\Users\\hp39wasi\\Dropbox\\sWorm\\SoilGrids\\1km"

  
  

sg.ftp <- "ftp://ftp.soilgrids.org/data/aggregated/1km/"
filenames = getURL(sg.ftp, ftp.use.epsv = TRUE, dirlistonly = TRUE)
filenames = strsplit(filenames, "\r*\n")[[1]]
filenames[1:5]


dl <- filenames[grep("PHIHOX", filenames)]
dl <- c(dl, filenames[grep("CLYPPT", filenames)])
dl <- c(dl, filenames[grep("SLTPPT", filenames)])
dl <- c(dl, filenames[grep("SNDPPT", filenames)])
dl <- c(dl, filenames[grep("CECSOL", filenames)])
dl <- c(dl, filenames[grep("ORCDRC", filenames)])
dl <- c(dl, filenames[grep("TAXNWRB_1", filenames)])

for(i in dl){
  try(download.file(file.path(sg.ftp, i), file.path(savedir, i)))
}
