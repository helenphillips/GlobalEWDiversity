loadMostRecent <- function(searchTerm, folderPath){
  files <- list.files(file.path(folderPath))
  files <- files[grep(searchTerm, files)]
  file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
  file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
  
  file_dates <- as.Date(file_dates)
  date <- max(file_dates, na.rm = TRUE)
  loadin <- files[grep(date, files)]

  return(loadin)
}

loadMostRecent_2 <- function(searchTerm, folderPath){
  files <- list.files(file.path(folderPath))
  files <- files[grep(searchTerm, files)]
  file_dates <- sapply(strsplit(files, "_"), "[", 3) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
  file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
  
  file_dates <- as.Date(file_dates)
  date <- max(file_dates, na.rm = TRUE)
  loadin <- files[grep(date, files)]
  
  return(loadin)
}
