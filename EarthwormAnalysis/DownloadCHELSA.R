library(RCurl)
library(stringr)

bioclim_int <- "https://www.wsl.ch/lud/chelsa/data/bioclim/integer"
savedir <- "C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\BioClim"

filenames = getURL(bioclim_int, ftp.use.epsv = FALSE, dirlistonly = TRUE)



html <- paste(readLines(bioclim_int), collapse="\n")
matched <- str_match_all(html, "<a href=\"(.*?)\"")
all <- matched[[1]][,2]
all <- all[grep("CHELSA", all)]

 for (filename in all) {
download.file(file.path(bioclim_int, filename), method="curl", file.path(savedir, filename), ssl.verifypeer = FALSE)
 }

####################################

prec <- "https://www.wsl.ch/lud/chelsa/data/climatologies/prec/"
savedir <- "C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\precip"

filenames = getURL(prec, ftp.use.epsv = FALSE, dirlistonly = TRUE)


html <- paste(readLines(prec), collapse="\n")
matched <- str_match_all(html, "<a href=\"(.*?)\"")
all <- matched[[1]][,2]
all <- all[grep("CHELSA", all)]

for (filename in all) {
  download.file(file.path(prec, filename), file.path(savedir, filename), mode = "w")
}

#####################################################
folders <- c("temp", "tmax", "tmin")

temp <- "https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/"
savedir <- "C:\\Users\\hp39wasi\\Dropbox\\sWorm\\CHELSAData\\temp"

for(f in folders){
  f_temp <- file.path(temp, f)
  f_savedir <- file.path(savedir, f)
  
  filenames = getURL(f_temp, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  
  
  html <- paste(readLines(f_temp), collapse="\n")
  matched <- str_match_all(html, "<a href=\"(.*?)\"")
  all <- matched[[1]][,2]
  all <- all[grep("CHELSA", all)]
  
  for (filename in all) {
    download.file(file.path(f_temp, filename), file.path(f_savedir, filename), mode = "w")
  }
  
}