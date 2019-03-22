########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}


data_out <- "Authorship"
#################################################
# 1. Load in models
#################################################

models <- "Models"


load(file.path(models, "richnessmodel.rds"))
load(file.path(models, "biomassmodel_full.rds"))
load(file.path(models, "abundancemodel_full.rds"))


#################################################
# 2. Loading in variables
#################################################


data_in <- "8_Data"


#################################################
# 3. Load in data
#################################################
files <- list.files(file.path(data_in))
files <- files[grep("sitesRichness_2", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

richness <- read.csv(file.path(data_in, loadin))

files <- list.files(file.path(data_in))
files <- files[grep("sitesBiomass_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

biomass <- read.csv(file.path(data_in, loadin))

files <- list.files(file.path(data_in))
files <- files[grep("sitesAbundance_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

abundance <- read.csv(file.path(data_in, loadin))

names(richness)[names(richness) == "scaleElevation"] <- "ScaleElevation"

alldat <- rbind(richness, biomass, abundance)

## unique sites
length(unique(alldat$ID))
## 7051
length(unique(alldat$Study_site))
## 7048

###############################################
## 4. Load in the bib data
################################################

bibdata <- "0_Data"

files <- list.files(file.path(bibdata))
files <- files[grep("Metadata_", files)]
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date
file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]

bib <- read.csv(file.path(bibdata, loadin))


files <- unique(alldat$file)

bib <- droplevels(bib[bib$file %in% files,])

length(unique(alldat$file)) == nrow(bib)
## 181

###################################################
## 5. Data Providers
###################################################

res <- data.frame(firstname = bib$DataProvider_FirstName, surname = bib$DataProvider_Surname, 
                  email = bib$DataProvider_Email)

# Get rid of rows with ALL NAS (I HAVE MESSED UP SOME DATA ENTRY)
res <- res[apply(res, 1, function(y) !all(is.na(y))),]
res$source <- "Data Provider"

naEmail <- res[is.na(res$email),]
Email <- res[!(is.na(res$email)),]
Email <- Email[order(Email$email),]


Email<-Email[!duplicated(Email$email),]

res <- rbind(naEmail, Email)
##
###################################################
## 5. DAta From Paper People
###################################################

res2 <- data.frame(firstname = NA, surname = bib$PaperContact_Surname, 
                  email = bib$PaperContact_Email)

naEmail <- res2[is.na(res2$email),]
naEmail <- naEmail[10,] ## The only one that has soem information

Email <- res2[!(is.na(res2$email)),]
Email <- Email[order(Email$email),]
Email<-Email[!duplicated(Email$email),]

res2 <- rbind(naEmail, Email)

res2$source <- "data from paper"
###################################################
## 5. Additional authors
###################################################


additional <- bib
additional <- additional[which(additional$Additional_Authors != "FirstName Surnam <emailaddress@somewhere>;"),]
additional <- additional[which(additional$Additional_Authors != "FirstName Surname <emailaddress@somewhere>;"),]

additional[,c(1:15)] <- NULL
additional[,c(2:ncol(additional))] <- NULL

# library(tidyr)
additional <- tidyr::separate(additional, 'Additional_Authors', paste("Author", 1:6, sep="_"), sep=";", extra="merge")

additional <- data.frame(authors = c(additional$Author_1, additional$Author_2,additional$Author_3,
                                     additional$Author_4,additional$Author_5,additional$Author_6))

additional <- as.data.frame(additional[complete.cases(additional),])

additional <- tidyr::separate(additional, 'additional[complete.cases(additional), ]', c('name', 'email'), sep="<", extra="merge")


### If name is blank we can remove the row
additional <- additional[which(additional$name != ""),]
additional<-additional[!duplicated(additional$name),]

### Remove the > at the end of the email
additional$email <- gsub(">", "", additional$email)
additional$email <- gsub("mailto:", "", additional$email)
additional$email <- gsub(" ", "", additional$email)

## Remove whitespace in name
additional$name <- trimws(additional$name, which = c("both"))


tmp <- data.frame(do.call(rbind, strsplit(additional$name, ' (?=[^ ]+$)', perl=TRUE)))

res3 <- cbind(additional, tmp)

res3 <- res3[,c('X1', 'X2', 'email')]
res3$source <- "additional author"
names(res3) <- names(res)



allResults <- rbind(res, res2, res3)

allResults$email <- gsub(">", "", allResults$email)
allResults$email <- gsub("<", "", allResults$email)



### Split by email address or not, so we can remove duplicates
allResults_noemail <- allResults[which(is.na(allResults$email)),]

allResults_email <- allResults[which(!(is.na(allResults$email))),]

#allResults_email <- allResults_email[order(allResults_email$email, allResults$firstname),]
allResults_email <- with(allResults_email, allResults_email[order(email,firstname),])

allResults_email<-allResults_email[!duplicated(allResults_email$email),]


allResults <- rbind(allResults_email, allResults_noemail)


##### SAVE

write.csv(allResults, file = file.path(data_out, "ListofDataProviders.csv"), row.names = FALSE)
## I then manually went through and highlighted those that are sWorm members