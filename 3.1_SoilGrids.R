library(XML)
library(RCurl)


j <- getURL("https://rest.soilgrids.org/query?lon=5.39&lat=51.57")


library(jsonlite)
j <- fromJSON("https://rest.soilgrids.org/query?lon=5.39&lat=51.57")


j1 <- fromJSON("https://rest.soilgrids.org/query?lon=5.39&lat=51.57&attributes=ORCDRC,CEC")

j2 <- fromJSON("https://rest.soilgrids.org/query?lon=5.39&lat=51.57&attributes=ORCDRC,CEC,CLYPPT,SNDPPT,SLTPPT,PHIHOX,TAXGWRB")
