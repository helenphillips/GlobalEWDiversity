library(XML)
library(RCurl)

sg.ftp <- "ftp://ftp.soilgrids.org/data/recent/"
filenames = getURL(sg.ftp, ftp.use.epsv = TRUE, dirlistonly = TRUE)
filenames = strsplit(filenames, "\r*\n")[[1]]
> filenames[1:5]
[1] "ACDWRB_M_ss_250m_ll.tif"     "ACDWRB_M_ss_250m_ll.tif.xml" "AWCh1_M_sl1_250m_ll.tif"
[4] "AWCh1_M_sl1_250m_ll.tif.xml" "AWCh1_M_sl2_250m_ll.tif"